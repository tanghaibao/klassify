#!/usr/bin/env python3
"""
Simulate PacBio reads with pbsim and generate CCS with pbccs (ccs).

Equivalent to:

  GENOME=$1
  PREFIX=`basename $GENOME .fa`

  pbsim --strategy wgs --method errhmm \
    --errhmm ~/code/pbsim3/data/ERRHMM-SEQUEL.model \
    --depth 50 --pass-num 10 --length-mean 15000 \
    --prefix $PREFIX --id-prefix $PREFIX --genome $GENOME

  ccs ${PREFIX}_0001.bam ${PREFIX}_0001.ccs.fastq.gz \
    --report-file ccs_report.txt --log-file ccs.log
"""

from __future__ import annotations
import argparse
import os
import sys

from pathlib import Path

from utils import check_tool, run

from jcvi.apps.base import cleanup, logger


def main():
    p = argparse.ArgumentParser(description="Run pbsim + ccs pipeline.")
    p.add_argument(
        "genome", type=Path, help="Reference genome FASTA (e.g., parents.genome.fa)"
    )
    p.add_argument(
        "--errhmm",
        default="~/code/pbsim3/data/ERRHMM-SEQUEL.model",
        help="ERRHMM model for pbsim (default: %(default)s)",
    )
    p.add_argument(
        "--depth", type=int, default=50, help="pbsim --depth (default: %(default)s)"
    )
    p.add_argument(
        "--pass-num",
        type=int,
        default=10,
        help="pbsim --pass-num (default: %(default)s)",
    )
    p.add_argument(
        "--length-mean",
        type=int,
        default=15000,
        help="pbsim --length-mean (default: %(default)s)",
    )
    p.add_argument(
        "--prefix",
        default=None,
        help="Output prefix (defaults to basename(GENOME) without .fa)",
    )
    p.add_argument(
        "--pbsim",
        default="pbsim",
        help="pbsim executable name/path (default: %(default)s)",
    )
    p.add_argument(
        "--ccs", default="ccs", help="pbccs executable name/path (default: %(default)s)"
    )
    args = p.parse_args()

    # Checks
    if not args.genome.exists():
        sys.exit(f"ERROR: Genome not found: {args.genome}")

    genome_dir = args.genome.parent.resolve()
    cwd = os.getcwd()
    os.chdir(genome_dir)

    check_tool(args.pbsim)
    check_tool(args.ccs)

    # Emulate: PREFIX=`basename $GENOME .fa`
    if args.prefix is not None:
        prefix = args.prefix
    else:
        name = args.genome.name
        prefix = (
            name[:-3] if name.endswith(".fa") else Path(name).stem
        )  # faithful + fallback

    errhmm_path = os.path.expanduser(args.errhmm)

    # Run pbsim
    run(
        [
            args.pbsim,
            "--strategy",
            "wgs",
            "--method",
            "errhmm",
            "--errhmm",
            errhmm_path,
            "--depth",
            str(args.depth),
            "--pass-num",
            str(args.pass_num),
            "--length-mean",
            str(args.length_mean),
            "--prefix",
            prefix,
            "--id-prefix",
            prefix,
            "--genome",
            str(args.genome),
        ]
    )

    # There are a few files that pbsim generates that we don't need
    cleanup(f"{prefix}_0001.maf.gz", f"{prefix}_0001.ref")

    # Run ccs on ${PREFIX}_0001.bam
    bam_in = f"{prefix}_0001.bam"
    ccs_out = f"{prefix}_0001.ccs.fastq.gz"
    run([args.ccs, bam_in, ccs_out])

    cleanup(f"{prefix}_0001.ccs.zmw_metrics.json.gz")

    logger.info("Done.")
    logger.info("  BAM:    %s", bam_in)
    logger.info("  CCS:    %s", ccs_out)
    os.chdir(cwd)


if __name__ == "__main__":
    main()
