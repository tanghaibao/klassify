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
    p.add_argument("--out", default="results", help="Directory to store the results.")
    p.add_argument(
        "-n", type=int, default=2, help="Number of genomes to use (2, 4, or 8)."
    )
    p.add_argument(
        "--seed", type=int, default=42, help="Random seed for reproducibility."
    )

    args = p.parse_args()
    genomes = get_genomes(args.ref)
    res = admix(genomes, args.n, seed=args.seed)
    print(res)


if __name__ == "__main__":
    main()
