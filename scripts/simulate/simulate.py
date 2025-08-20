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
- Run klassify pipeline with the generated resources
- Move results to the results directory
"""

import argparse
import os
import os.path as op
import random

from dataclasses import dataclass
from itertools import batched
from random import sample
from typing import List

from jcvi.apps.base import logger, mkdir
from jcvi.formats.base import FileMerger


@dataclass
class SimulatedConfig:
    parents: List[str]
    gametes: List[str]
    seed: int


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


def prepare_input(
    config: SimulatedConfig, ref_dir: str, mosaics_dir: str, out_dir: str
):
    """
    Prepare the input files for the klassify pipeline.
    """
    mkdir(out_dir)
    parents_genomes = op.join(out_dir, "parents.genomes.fa")
    parent_reads = op.join(out_dir, "parent_reads.fq.gz")
    f1_reads = op.join(out_dir, "f1_reads.fq.gz")
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
    p.add_argument("--out", default="out", help="Directory to store the results.")
    p.add_argument(
        "-n", type=int, default=2, help="Number of genomes to use (2, 4, or 8)."
    )
    p.add_argument(
        "--seed", type=int, default=42, help="Random seed for reproducibility."
    )

    args = p.parse_args()
    genomes = get_genomes(args.ref)
    res = admix(genomes, args.n, seed=args.seed)
    logger.info(res)

    # Prepare input files for the klassify pipeline
    prepare_input(res, args.ref, args.mosaics, args.out)
    logger.info("Input files prepared in `%s`", args.out)


if __name__ == "__main__":
    main()
