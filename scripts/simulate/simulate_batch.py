#!/usr/bin/env python3
"""
Run simulate.py on:
  - ref
  - mosaics
in parallel.

Usage:
  python build_reads_batch.py --jobs 16
"""

import argparse
import os
import sys

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Tuple

from utils import run

from jcvi.apps.base import logger


def run_one(script: Path, n_seed: Tuple[int, int]) -> Tuple[Tuple[int, int], int]:
    """
    Run build_reads.py on a single FASTA path, returning (path, returncode).
    """
    n, seed = n_seed
    cmd = [sys.executable, str(script), "-n", str(n), "--seed", str(seed), "--cleanup"]
    # Inherit stdout/stderr so you can see script output live per task
    proc = run(cmd)
    return n_seed, proc.returncode


def main():
    parser = argparse.ArgumentParser(description="Parallel runner for build_reads.py")
    parser.add_argument(
        "--jobs",
        "-j",
        type=int,
        default=24,
        help="Number of parallel workers (default: %(default)d)",
    )
    parser.add_argument(
        "--script",
        type=str,
        default="~/code/klassify/scripts/simulate/simulate.py",
        help="Path to simulate.py (default: %(default)s)",
    )

    args = parser.parse_args()

    script_path = Path(os.path.expanduser(args.script))
    if not script_path.exists():
        logger.error("script not found: %s", script_path)
        sys.exit(2)

    logger.info("[INFO] Using script: %s", script_path)
    logger.info("[INFO] Jobs: %d", args.jobs)
    failures = []

    # ThreadPool is ideal for launching many subprocesses
    Ns = [2, 4, 8]
    seeds = range(100)
    total_jobs = len(Ns) * len(seeds)
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        n_seeds, futures = {
            ex.submit(run_one, script_path, (n, seed)): (n, seed)
            for n in Ns
            for seed in seeds
        }
        for fut in as_completed(futures):
            rc = fut.result()
            if rc == 0:
                logger.info("[OK] %s", n_seeds)
            else:
                logger.info("[FAIL:%d] %s", rc, n_seeds)
                failures.append((n_seeds, rc))

    if failures:
        logger.info(
            "[SUMMARY] %d failures out of %d tasks:",
            len(failures),
            total_jobs,
        )
        for n_seeds, rc in failures:
            logger.info("  - %s (exit %d)", n_seeds, rc)
        sys.exit(1)
    else:
        logger.info("[SUMMARY] All %d tasks completed successfully.", total_jobs)
        sys.exit(0)


if __name__ == "__main__":
    main()
