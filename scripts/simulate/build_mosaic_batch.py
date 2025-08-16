#!/usr/bin/env python3
"""
Run build_mosaic.py for all ordered pairs in ref/*.fasta, each in its own dir.

Example
-------
python build_mosaic_batch.py \
  --ref-dir ref \
  --build build_mosaic.py \
  --out-root . \
  --jobs 8 \
  --extra "--min-distance 0 --preset map-hifi"

This creates A_B/, A_C/, ... and runs:
  python build_mosaic.py <A.fa> <B.fa> [--extra ...]   (cwd = A_B/)
"""

from __future__ import annotations
import argparse
import os
import sys

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from build_mosaic import build_mosaic

from jcvi.apps.base import logger, mkdir


FA_EXTS = (".fa", ".fasta", ".fna")


def find_fastas(ref_dir: Path) -> list[Path]:
    """
    Find all FASTA files in ref_dir.
    """
    files = []
    for ext in FA_EXTS:
        files.extend(sorted(ref_dir.glob(f"*{ext}")))
    if not files:
        raise SystemExit(f"No FASTA files found in: {ref_dir} (looked for {FA_EXTS})")
    return files


def stem(p: Path) -> str:
    """
    Remove first matching extension from FA_EXTS
    """
    # remove first matching extension from FA_EXTS, else use .stem
    for ext in FA_EXTS:
        if p.name.endswith(ext):
            return p.name[: -len(ext)]
    return p.stem


def ordered_pairs(items: list[Path]) -> list[tuple[Path, Path]]:
    """
    Return all ordered pairs of items.
    """
    return [(a, b) for a in items for b in items if a != b]


def run_one_pair(
    fa_a: Path,
    fa_b: Path,
    out_root: Path,
    dry_run: bool,
    min_distance: int,
    n: int,
    seed: int,
) -> tuple[str, int]:
    """
    Run build_mosaic.py for a single pair of FASTA files.
    """
    A, B = stem(fa_a), stem(fa_b)
    pair = f"{A}_{B}"
    pair_dir = out_root / pair
    mkdir(pair_dir)

    cmd = [str(fa_a), str(fa_b), str(pair_dir / pair), None, 1, min_distance, n, seed]

    if dry_run:
        return (" ".join(str(x) for x in cmd) + f"   (cwd={pair_dir})", 0)

    cwd = os.getcwd()
    os.chdir(pair_dir)
    build_mosaic(*cmd)
    os.chdir(cwd)
    return (f"[OK] {pair_dir}", 0)


def main():
    ap = argparse.ArgumentParser(
        description="Batch runner for build_mosaic.py over ordered pairs."
    )
    ap.add_argument(
        "ref_dir",
        type=Path,
        help="Directory containing *.fasta (8 files).",
    )
    ap.add_argument(
        "--out-root",
        type=Path,
        default=Path("mosaics"),
        help="Where to create A_B/ dirs (default: %(default)s)",
    )
    ap.add_argument(
        "--jobs", type=int, default=8, help="Parallel jobs (default: %(default)s)"
    )
    ap.add_argument("--dry-run", action="store_true", help="Print commands only")
    ap.add_argument(
        "--unordered",
        action="store_true",
        help="Use unordered pairs (nC2) instead of ordered (n*(n-1))",
    )
    ap.add_argument(
        "--min-distance",
        type=int,
        default=1_000_000,
        help="Minimum spacing between A breakpoints (default: %(default)s bp)",
    )
    ap.add_argument(
        "--n", type=int, default=4, help="Number of breakpoints to simulate"
    )
    ap.add_argument("--seed", type=int, default=42)

    args = ap.parse_args()
    ref_dir = args.ref_dir.resolve()
    out_root = args.out_root.resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    fastas = find_fastas(ref_dir)

    # Determine pairs
    if args.unordered:
        # nC2
        pairs = []
        for i in range(len(fastas)):
            for j in range(i + 1, len(fastas)):
                pairs.append((fastas[i], fastas[j]))
    else:
        # ordered n*(n-1)  → 8*7 = 56
        pairs = ordered_pairs(fastas)

    logger.info("Found %d FASTA(s) in `%s`", len(fastas), ref_dir)
    logger.info(
        "Planning %d %s pair runs.",
        len(pairs),
        "unordered" if args.unordered else "ordered",
    )
    if args.dry_run:
        logger.info("DRY RUN (no commands will execute)\n")

    # Parallel execution
    rc_total = 0
    seed = args.seed
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futs = []
        for fa_a, fa_b in pairs:
            fut = ex.submit(
                run_one_pair,
                fa_a,
                fa_b,
                out_root,
                args.dry_run,
                args.min_distance,
                args.n,
                seed,
            )
            futs.append(fut)
            seed += 1  # increment seed for each pair

        for fut in as_completed(futs):
            msg, rc = fut.result()
            logger.info(msg)
            rc_total |= rc

    if rc_total != 0:
        sys.exit(rc_total)


if __name__ == "__main__":
    main()
