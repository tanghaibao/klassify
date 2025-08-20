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

def main():
    pass


if __name__ == "__main__":
    main()