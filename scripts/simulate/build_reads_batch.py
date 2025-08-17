#!/usr/bin/env python3
"""
Run build_reads.py on:
  - ref/*.fa
  - mosaics/*/*.mosaic.fa
in parallel.

Usage:
  python build_reads_batch.py --jobs 8
"""
from __future__ import annotations
import argparse
import os
import sys

from pathlib import Path
from glob import glob
from concurrent.futures import ThreadPoolExecutor, as_completed

from utils import run
from jcvi.apps.base import logger


def gather_inputs(patterns: list[str]) -> list[str]:
    """
    Gather input files from the given patterns.
    """
    files: list[str] = []
    for pat in patterns:
        files.extend(glob(os.path.expanduser(pat)))
    # absolute + dedup + sort for stable ordering
    files = sorted(set(os.path.abspath(p) for p in files))
    return files


def run_one(script: Path, fasta_path: str) -> tuple[str, int]:
    """
    Run build_reads.py on a single FASTA path, returning (path, returncode).
    """
    cmd = [sys.executable, str(script), fasta_path]
    # Inherit stdout/stderr so you can see script output live per task
    proc = run(cmd)
    return (fasta_path, proc.returncode)


def main():
    parser = argparse.ArgumentParser(description="Parallel runner for build_reads.py")
    parser.add_argument(
        "--jobs",
        "-j",
        type=int,
        default=16,
        help="Number of parallel workers (default: %(default)d)",
    )
    parser.add_argument(
        "--script",
        type=str,
        default="~/code/klassify/scripts/simulate/build_reads.py",
        help="Path to build_reads.py (default: %(default)s)",
    )
    parser.add_argument(
        "--pattern",
        action="append",
        dest="patterns",
        help="Glob pattern(s) to search (can be given multiple times). "
        "If omitted, uses: 'ref/*.fa' and 'mosaics/*/*.mosaic.fa'.",
    )
    args = parser.parse_args()

    script_path = Path(os.path.expanduser(args.script))
    if not script_path.exists():
        logger.error("[ERROR] script not found: %s", script_path)
        sys.exit(2)

    patterns = args.patterns or ["ref/*.fa", "mosaics/*/*.mosaic.fa"]
    inputs = gather_inputs(patterns)

    if not inputs:
        logger.warning("[WARN] No input files matched the given patterns.")
        sys.exit(0)

    logger.info("[INFO] Using script: %s", script_path)
    logger.info("[INFO] Jobs: %d", args.jobs)
    logger.info("[INFO] Files to process: %d", len(inputs))

    failures = []
    # ThreadPool is ideal for launching many subprocesses
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futures = {ex.submit(run_one, script_path, f): f for f in inputs}
        for fut in as_completed(futures):
            path, rc = fut.result()
            if rc == 0:
                logger.info("[OK] %s", path)
            else:
                logger.info("[FAIL:%d] %s", rc, path)
                failures.append((path, rc))

    if failures:
        logger.info(
            "[SUMMARY] %d failures out of %d tasks:", len(failures), len(inputs)
        )
        for path, rc in failures:
            logger.info("  - %s (exit %d)", path, rc)
        sys.exit(1)
    else:
        logger.info("[SUMMARY] All %d tasks completed successfully.", len(inputs))
        sys.exit(0)


if __name__ == "__main__":
    main()
