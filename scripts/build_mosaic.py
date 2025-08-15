#!/usr/bin/env python3
"""
Breakpoint-first mosaic using a precomputed UCSC chain (minimap2 + transanno already run),
with optional on-the-fly generation if no chain is provided.

Pipeline:
1) (Optional) If --chain not provided:
     minimap2 (A vs B) -> PAF
     transanno minimap2-to-chain -> CHAIN
2) Load CHAIN with the 'liftover' Python library to get an A->B converter.
3) Simulate N breakpoints on A with a minimum spacing.
4) Lift each A breakpoint to B (keep only those mapping to B's main chrom on '+'
   strand; require strictly increasing on B).
5) Write paired A/B breakpoints to TSV.
6) Build mosaic by alternating A and B between consecutive breakpoints, starting on A.

Requirements:
  - Python: liftover, pyfaidx, numpy
  - If --chain is NOT supplied: minimap2 and transanno must be in PATH.

Example:
  python mosaic_liftover.py \
    --faA A.fa --faB B.fa \
    --out-prefix results/run1 \
    --chain precomputed.chain \
    --n 4 --min-distance 1000000 --seed 13
"""

import argparse
import os
import shutil
import subprocess

from typing import List, Tuple

import numpy as np

from liftover import ChainFile
from pyfaidx import Fasta


def read_first_record(fa_path: str) -> Tuple[str, str]:
    """
    Read the first record from a FASTA file.
    """
    fa = Fasta(fa_path)
    if not fa.keys():
        raise ValueError(f"No records in {fa_path}")
    name = list(fa.keys())[0]
    return name, str(fa[name][:].seq).upper()


def write_fasta(path: str, header: str, seq: str, width: int = 60):
    """
    Write a FASTA file.
    """
    with open(path, "w", encoding="utf-8") as fw:
        fw.write(f">{header}\n")
        for i in range(0, len(seq), width):
            fw.write(seq[i : i + width] + "\n")


def ensure_tool(name: str):
    """
    Ensure a tool is in the PATH.
    """
    if shutil.which(name) is None:
        raise RuntimeError(f"Required tool '{name}' not found in PATH")


def run_minimap2_paf(faA: str, faB: str, threads: int, out_paf: str):
    """
    Produce a PAF mapping A (query) -> B (target) with CIGAR/CS; suitable for
    transanno.
    """
    cmd = ["minimap2", "-cx", "asm5", "--cs", "-t", str(threads), faB, faA]
    with open(out_paf, "w", encoding="utf-8") as out:
        subprocess.run(cmd, check=True, stdout=out)


def paf_to_chain_with_transanno(paf: str, out_chain: str):
    """
    Convert PAF to UCSC chain (A -> B)
    """
    cmd = ["transanno", "minimap2chain", paf, "--output", out_chain]
    subprocess.run(cmd, check=True)


def sample_breakpoints(L: int, n: int, dmin: int, seed: int) -> List[int]:
    """
    Sample n positions in [1, L-1) (avoid exact ends), with >= dmin spacing.
    Greedy retries; may return fewer if space is tight.
    """
    if L <= 2 or n <= 0:
        return []
    rng = np.random.default_rng(seed)
    picks: List[int] = []
    tries = 0
    max_tries = 200000
    low, high = 1, L - 1
    while len(picks) < n and tries < max_tries:
        x = int(rng.integers(low, high))
        if all(abs(x - y) >= dmin for y in picks):
            picks.append(x)
        tries += 1
    picks.sort()
    return picks


def lift_pos(converter, chromA: str, posA_0based: int):
    """
    Use liftover to map a 0-based A position to B.
    Returns list of hits (chrom, pos, strand); we will filter outside.
    """
    return converter[chromA][posA_0based]


def choose_mapped_breaks(
    converter, chromA: str, chromB_target: str, L: int, n: int, dmin: int, seed: int
) -> Tuple[List[int], List[int]]:
    """
    Find n A breakpoints that:
      - are >= dmin apart on A,
      - each lifts to B (same chrom as chromB_target),
      - all on '+' strand (to keep B slicing simple),
      - B positions strictly increase with A.
    Will resample multiple times; raises if cannot satisfy.
    """
    rng = np.random.default_rng(seed)
    attempts = 200
    for _ in range(attempts):
        a_breaks = sample_breakpoints(L, n, dmin, int(rng.integers(0, 2**31 - 1)))
        if len(a_breaks) < n:
            continue
        b_hits: List[int] = []
        ok = True
        prev_b = -1
        for a in a_breaks:
            candidates = lift_pos(converter, chromA, a)
            # Filter candidates: same target chrom, plus strand
            candidates = [
                (c, p, s) for (c, p, s) in candidates if c == chromB_target and s == "+"
            ]
            if not candidates:
                ok = False
                break
            # Choose the first (liftover may return multiple blocks)
            _, bpos, _ = candidates[0]
            bpos = int(bpos)
            if bpos <= prev_b:
                ok = False
                break
            b_hits.append(bpos)
            prev_b = bpos
        if ok:
            return a_breaks, b_hits
    raise RuntimeError(
        "Could not find enough mapped breakpoints; try lowering --min-distance or using a different seed."
    )


def build_mosaic_alternate(
    seqA: str, seqB: str, a_breaks: List[int], b_breaks: List[int]
) -> str:
    """
    Start on A. Between consecutive A breakpoints, alternate A and B segments:
      segments: [0..a1] (A), [b1..b2] (B), [a2..a3] (A), ...
    Tail: append the remainder of whichever donor you end on.
    """
    assert len(a_breaks) == len(b_breaks) and len(a_breaks) > 0
    L = len(seqA)
    a_pos = [0] + a_breaks + [L]
    pieces: List[str] = []
    donor_A = True  # first interior segment is A
    for i in range(len(a_pos) - 1):
        aL, aR = a_pos[i], a_pos[i + 1]
        if donor_A:
            pieces.append(seqA[aL:aR])
        else:
            # B segment uses corresponding b[i-1]..b[i]
            bi = i - 1
            bL, bR = b_breaks[bi], b_breaks[bi + 1]
            pieces.append(seqB[bL:bR] if bL <= bR else "")
        donor_A = not donor_A
    return "".join(pieces)


def main():
    ap = argparse.ArgumentParser(
        description="Breakpoint-first mosaic using a precomputed chain (or generate one if absent)."
    )
    ap.add_argument(
        "--faA", required=True, help="FASTA for genome A (source of breakpoints)"
    )
    ap.add_argument(
        "--faB", required=True, help="FASTA for genome B (target of liftover)"
    )
    ap.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix (creates .tsv and .fa; may also create .paf/.chain if --chain not provided)",
    )
    ap.add_argument(
        "--chain",
        help="Precomputed UCSC chain file mapping A->B. If provided, minimap2/transanno are skipped.",
    )
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument(
        "--min-distance",
        type=int,
        default=1_000_000,
        help="Minimum spacing between A breakpoints (bp)",
    )
    ap.add_argument(
        "--n", type=int, default=4, help="Number of breakpoints to simulate"
    )
    ap.add_argument("--seed", type=int, default=13)
    args = ap.parse_args()

    chromA, seqA = read_first_record(args.faA)
    chromB, seqB = read_first_record(args.faB)

    out_tsv = args.out_prefix + ".breakpoints.tsv"
    out_fa = args.out_prefix + ".mosaic.fa"

    # Get or make chain
    if args.chain:
        chain_path = args.chain
        if not os.path.exists(chain_path):
            raise FileNotFoundError(f"--chain file not found: {chain_path}")
    else:
        ensure_tool("minimap2")
        ensure_tool("transanno")
        out_paf = args.out_prefix + ".paf"
        chain_path = args.out_prefix + ".chain"
        run_minimap2_paf(args.faA, args.faB, args.threads, out_paf)
        paf_to_chain_with_transanno(out_paf, chain_path)

    # Load liftover converter
    converter = ChainFile(chain_path, one_based=False)

    # Simulate and map breakpoints
    a_breaks, b_breaks = choose_mapped_breaks(
        converter=converter,
        chromA=chromA,
        chromB_target=chromB,
        L=len(seqA),
        n=int(args.n),
        dmin=int(args.min_distance),
        seed=int(args.seed),
    )

    # Write TSV of A/B breakpoints (0-based)
    os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
    with open(out_tsv, "w", encoding="utf-8") as w:
        w.write("idx\tA_chrom\tA_bp_0based\tB_chrom\tB_bp_0based\tstrand\n")
        for i, (a, b) in enumerate(zip(a_breaks, b_breaks)):
            w.write(f"{i}\t{chromA}\t{a}\t{chromB}\t{b}\t+\n")

    # Build mosaic (alternate A/B between breakpoints)
    mosaic = build_mosaic_alternate(seqA, seqB, a_breaks, b_breaks)
    write_fasta(out_fa, header=os.path.basename(args.out_prefix), seq=mosaic)


if __name__ == "__main__":
    main()
