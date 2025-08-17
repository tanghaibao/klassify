#!/usr/bin/env python3

"""
Extract chr1 from each input FASTA(.gz) and write to ref/<prefix>.fasta with a new ID.

Defaults:
- Input glob: dir/*.fasta.gz
- Output dir: ref/
- New header: >{prefix}.chr1   (prefix = basename without .fasta[.gz])

Matching:
- By default matches common chr1 spellings:
  - "chr1", "Chr1", "chr01", "Chr01", "1" (as a whole token),
  - "chromosome 1", "chromosome_1", etc.
- Override with --pattern if needed (regex, applied to the first header token).

Examples
--------
# Simple run, 8 files in dir/*.fasta.gz → 8 outputs in ref/
python extract_chr1_to_ref.py

# Custom input glob and more threads
python extract_chr1_to_ref.py --glob "genomes/*.fa.gz" --jobs 4

# If your headers are unusual, specify a custom regex:
python extract_chr1_to_ref.py --pattern "(?i)^NC_003070(\.|$)"   # TAIR10 Chr1 RefSeq
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable, Tuple

from jcvi.apps.base import logger

FA_EXTS = (".fa", ".fasta", ".fna")

DEFAULT_PATTERN = (
    r"(?ix)"  # ignorecase, verbose
    r"^("  # start of token
    r"chr[_\-\s]*0*1"  # chr1, chr01, chr-1, chr_1
    r"|chromosome[_\-\s]*0*1"  # chromosome 1
    r"|0*1"  # plain '1' (but as a whole token)
    r")$"
)


def parse_args():
    ap = argparse.ArgumentParser(
        description="Extract chr1 → ref/<prefix>.fasta with renamed header."
    )
    ap.add_argument(
        "--glob",
        default="version2.5.2019-10-09/*.fasta.gz",
        help="Input glob for FASTA files (default: %(default)s)",
    )
    ap.add_argument(
        "--outdir",
        type=Path,
        default=Path("ref"),
        help="Output directory (default: %(default)s)",
    )
    ap.add_argument(
        "--pattern",
        default=DEFAULT_PATTERN,
        help="Regex applied to the first header token to detect chr1 (default matches common chr1 spellings)",
    )
    ap.add_argument(
        "--id-sep",
        default="_",
        help="Separator between prefix and chr1 in new header (default: '%(default)s')",
    )
    ap.add_argument(
        "--jobs", type=int, default=8, help="Parallel workers (default: %(default)d)"
    )
    ap.add_argument("--force", action="store_true", help="Overwrite existing outputs")
    ap.add_argument(
        "--dry-run", action="store_true", help="Print planned actions without writing"
    )
    return ap.parse_args()


def prefix_until_dot_or_underscore(s: str) -> str:
    """
    Get the prefix from a string, up to the first dot or underscore.
    """
    return s.split(".")[0].split("_")[0]


def guess_prefix(p: Path) -> str:
    """
    Guess the prefix from the file name.
    """
    name = p.name
    if name.endswith(".gz"):
        name = name[:-3]
    for ext in FA_EXTS:
        if name.endswith(ext):
            return prefix_until_dot_or_underscore(name[: -len(ext)])
    s = Path(name).stem
    return prefix_until_dot_or_underscore(s)


def nopen(path: Path):
    """
    Open a file, automatically detecting .gz suffix.
    """
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "rt")


def iter_fasta_records(path: Path) -> Iterable[Tuple[str, Iterable[str]]]:
    """
    Yield (header_line_without_>, seq_lines_iterable) for each FASTA record.
    seq_lines_iterable yields raw sequence lines (no newlines).
    """
    handle = nopen(path)
    header = None
    buf = []
    try:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq_iter = (s for s in buf)
                    yield header, seq_iter
                header = line[1:].strip()
                buf = []
            else:
                buf.append(line)
        if header is not None:
            seq_iter = (s for s in buf)
            yield header, seq_iter
    finally:
        handle.close()


def first_token(header: str) -> str:
    """
    Get the first token from a header.
    """
    return header.split()[0]


def find_chr1_token(headers: Iterable[str], pattern: re.Pattern) -> str | None:
    """
    Find the first header that matches the pattern.
    """
    for h in headers:
        tok = first_token(h)
        if pattern.search(tok):
            return h
    return None


def extract_chr1(in_path: Path, out_path: Path, new_id: str, pat: re.Pattern) -> int:
    """
    Stream-extract the first record whose first token matches pat.
    Write as >new_id and the sequence lines as-is.
    Returns the number of bases written (0 if not found).
    """
    bases = 0
    target_header = None
    # First pass to find matching header quickly (no need, but gives clearer error)
    headers = []
    for h, _ in iter_fasta_records(in_path):
        headers.append(h)
        if pat.search(first_token(h)):
            target_header = h
            break

    if target_header is None:
        # No match; help user by listing a few tokens
        tokens = [first_token(h) for h in headers[:10]]
        raise RuntimeError(
            f"No chr1-like header found in {in_path}.\n"
            f"Seen tokens (first 10): {tokens}\n"
            "Consider adjusting --pattern."
        )

    # Second pass to write it
    with open(out_path, "wt") as out:
        out.write(f">{new_id}\n")
        for h, seq_iter in iter_fasta_records(in_path):
            if first_token(h) == first_token(target_header):
                # Write lines as in source (no rewrap)
                for s in seq_iter:
                    out.write(s + "\n")
                    bases += len(s)
                break
    return bases


def process_one(
    in_path: Path,
    outdir: Path,
    pat: re.Pattern,
    id_sep: str,
    force: bool,
    dry_run: bool,
) -> str:
    """
    Process a single FASTA file.
    """
    prefix = guess_prefix(in_path)
    out_path = outdir / f"{prefix}.fa"
    new_id = f"{prefix}{id_sep}chr1"

    if out_path.exists() and not force:
        return f"[SKIP] {out_path} exists"

    if dry_run:
        return f"[DRY] {in_path} → {out_path}   header: >{new_id}"

    outdir.mkdir(parents=True, exist_ok=True)
    n = extract_chr1(in_path, out_path, new_id, pat)
    return f"[OK] {in_path.name} → {out_path.name} ({n:,} bp)"


def main():
    args = parse_args()
    files = sorted(Path().glob(args.glob))
    if not files:
        sys.exit(f"No files matched: {args.glob}")

    pat = re.compile(args.pattern)

    logger.info("Found %d file(s). Writing to `%s/`", len(files), args.outdir)
    if args.dry_run:
        logger.info("DRY RUN (no files will be written)\n")

    results = []
    rc = 0
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futs = {
            ex.submit(
                process_one, f, args.outdir, pat, args.id_sep, args.force, args.dry_run
            ): f
            for f in files
        }
        for fut in as_completed(futs):
            try:
                msg = fut.result()
                logger.info(msg)
                results.append(msg)
            except Exception as e:
                rc = 1
                logger.error(f"[FAIL] {futs[fut].name}: %s", e)

    sys.exit(rc)


if __name__ == "__main__":
    main()
