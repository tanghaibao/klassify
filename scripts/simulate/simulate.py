#!/usr/bin/env python3
"""
This script kicks off the breakpoint runs given a set of artificial mossaic genomes.

Procedure:
- Find all the 8 genomes in the ref directory
- Choose N genomes at random, diploid (N=2), tetraploid (N=4), octoploid (N=8)
- Randomly select pairs of genomes to get mosaics
- Build the following resources:
    - parents.genomes.fa: concatenated N genomes
    - parent_reads.fq.gz: concatenated reads from the N genomes
    - f1_reads.fq.gz: concatenated reads from the mosaic genomes
    - breakpoints.tsv: concatenated breakpoints from the mosaic genomes
- Run klassify pipeline with the generated resources
- Move results to the results directory
"""

import argparse
import os
import os.path as op
import shutil

import pandas as pd
import random

from dataclasses import dataclass
from itertools import batched
from pathlib import Path
from random import sample
from typing import List

from jcvi.apps.base import logger, mkdir, sh
from jcvi.formats.base import FileMerger

DEFAULT_SCRIPT_PATH = Path(__file__).parent.parent
PARENTS_GENOMES = "parents.genomes.fa"
PARENT_READS = "parent_reads.fq.gz"
F1_READS = "f1_reads.fq.gz"
REGIONS_OUTPUT = "f1_classify.paired.regions"


@dataclass
class SimulatedConfig:
    parents: List[str]
    gametes: List[str]
    seed: int

    @property
    def name(self):
        """
        Generate a name for the simulation based on the parents and seed.
        """
        return f"simulated_{'+'.join(self.gametes)}_seed{self.seed}"


def get_genomes(ref_dir: str):
    """
    Get a list of genome files from the reference directory.
    """
    genomes = [
        op.basename(f).split(".", 1)[0]
        for f in os.listdir(ref_dir)
        if f.endswith(".fa")
    ]
    if len(genomes) < 8:
        raise ValueError("Not enough genomes found in the reference directory.")
    return genomes


def admix(genomes: List[str], n: int, seed: int) -> SimulatedConfig:
    """
    Create mosaic genomes by randomly selecting pairs of genomes.
    """
    random.seed(seed)
    parents = sample(genomes, n)
    gametes = ["_".join(x) for x in batched(parents, 2)]
    return SimulatedConfig(
        parents,
        gametes,
        seed,
    )


def prepare_input(config: SimulatedConfig, ref_dir: str, mosaics_dir: str):
    """
    Prepare the input files for the klassify pipeline.
    """
    out_dir = config.name
    mkdir(out_dir)
    parents_genomes = op.join(out_dir, PARENTS_GENOMES)
    parent_reads = op.join(out_dir, PARENT_READS)
    f1_reads = op.join(out_dir, F1_READS)
    parents = [op.join(ref_dir, f"{p}.fa") for p in config.parents]
    # Concatenate parent genomes
    FileMerger(parents, parents_genomes).merge()
    logger.info("Parents genomes written to `%s`", parents_genomes)
    # Concatenate parent reads
    parent_reads_list = [f"{ref_dir}/{p}_0001.ccs.fastq.gz" for p in config.parents]
    FileMerger(parent_reads_list, parent_reads).merge()
    logger.info("Parent reads written to `%s`", parent_reads)
    # Concatenate mosaic reads
    f1_reads_list = [
        f"{mosaics_dir}/{gamete}/{gamete}.mosaic_0001.ccs.fastq.gz"
        for gamete in config.gametes
    ]
    FileMerger(f1_reads_list, f1_reads).merge()
    logger.info("F1 reads written to `%s`", f1_reads)


def prepare_output(config: SimulatedConfig, mosaics_dir: str, results_dir: str):
    """
    Prepare the output files for the klassify pipeline.
    """
    out_dir = config.name
    mkdir(results_dir)
    # Concatenate true breakpoints
    breakpoints = []
    for gamete in config.gametes:
        bp_file = f"{mosaics_dir}/{gamete}/{gamete}.breakpoints.tsv"
        if op.exists(bp_file):
            breakpoints.append(pd.read_csv(bp_file, sep="\t"))
    df = pd.concat(breakpoints, ignore_index=True)
    df["A_breakpoint"] = df["A_chrom"] + ":" + df["A_bp_0based"].astype(str)
    df["B_breakpoint"] = df["B_chrom"] + ":" + df["B_bp_0based"].astype(str)
    df["Run"] = config.name
    df["Source"] = "true"
    columns = ["Run", "Source", "A_breakpoint", "B_breakpoint"]
    df = df[columns]

    # Now get the computed breakpoints
    regions_output = op.join(out_dir, REGIONS_OUTPUT)
    if not op.exists(regions_output):
        raise FileNotFoundError(f"Regions output file not found: {regions_output}")
    rows = [x.strip() for x in open(regions_output)]
    computed_df = []
    for a, b in batched(rows, 2):
        a_genome, b_genome = a.split("_", 1)[0], b.split("_", 1)[0]
        if f"{a_genome}_{b_genome}" not in config.gametes:
            a, b = b, a
        computed_df.append((a, b))
    computed_df = pd.DataFrame(computed_df, columns=["A_breakpoint", "B_breakpoint"])
    computed_df["Run"] = config.name
    computed_df["Source"] = "computed"
    computed_df = computed_df[columns]
    df = pd.concat([df, computed_df], ignore_index=True)

    breakpoints_output = op.join(results_dir, f"{config.name}.breakpoints.tsv")
    df.to_csv(breakpoints_output, sep="\t", index=False)
    logger.info("True and computed breakpoint written to `%s`", breakpoints_output)


def main():
    p = argparse.ArgumentParser(
        description="Simulate mosaic genomes and run the klassify pipeline."
    )
    p.add_argument(
        "--ref", default="ref", help="Directory containing reference genomes."
    )
    p.add_argument(
        "--mosaics", default="mosaics", help="Directory containing reference genomes."
    )
    p.add_argument(
        "--results-dir", default="results", help="Directory to store the results."
    )
    p.add_argument(
        "-n", type=int, default=2, help="Number of genomes to use (2, 4, or 8)."
    )
    p.add_argument(
        "--seed", type=int, default=42, help="Random seed for reproducibility."
    )
    p.add_argument(
        "--cleanup",
        default=False,
        action="store_true",
        help="Remove run dir after processing.",
    )
    p.add_argument(
        "--scripts-dir",
        type=Path,
        default=DEFAULT_SCRIPT_PATH,
        help="Directory containing helper scripts (pipeline.py).",
    )

    args = p.parse_args()
    genomes = get_genomes(args.ref)
    config = admix(genomes, args.n, seed=args.seed)
    logger.info(config)

    # Prepare input files for the klassify pipeline
    prepare_input(config, args.ref, args.mosaics)
    out_dir = config.name
    logger.info("Input files prepared in `%s`", out_dir)

    # Run klassify pipeline
    cwd = Path.cwd()
    os.chdir(cwd / out_dir)
    klassify_script = args.scripts_dir / "pipeline.py"
    if not klassify_script.exists():
        raise FileNotFoundError(f"Pipeline script not found: {klassify_script}")
    cmd = f"python {klassify_script} {F1_READS} {PARENT_READS} {PARENTS_GENOMES}"
    if not op.exists(REGIONS_OUTPUT):
        sh(cmd)
    os.chdir(cwd)
    logger.info("Klassify pipeline completed in `%s`.", out_dir)

    # Write results to the output directory
    prepare_output(
        config,
        args.mosaics,
        args.results_dir,
    )

    # Clean up the run directory if requested
    if args.cleanup:
        logger.info("Cleaning up run directory `%s`.", out_dir)
        shutil.rmtree(out_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
