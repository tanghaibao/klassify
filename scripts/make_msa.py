import sys

from typing import List, Set

import pandas as pd
import pysam

from Bio import pairwise2

from jcvi.apps.base import logger, sh


def write_aligned_fasta(region: str, ban: Set[str], aligned_fasta: str):
    """
    Write multiple sequence FASTA around breakpoints.
    """
    ra, rb = region.split()
    PAD = 500
    fasta = pysam.FastaFile(
        "/Users/bao/projects/female-restitution/ref/SoSs_ref2024.fa"
    )
    a, astart_aend = ra.split(":")
    b, bstart_bend = rb.split(":")
    _, aend = [int(x) for x in astart_aend.split("-")]
    bstart, _ = [int(x) for x in bstart_bend.split("-")]
    seqa = fasta.fetch(a, aend - PAD, aend + PAD)
    seqb = fasta.fetch(b, bstart - PAD, bstart + PAD)
    seqs = [
        (len(seqa), f"{a}:{aend - PAD}:{aend + PAD}", seqa),
        (len(seqb), f"{b}:{bstart - PAD}:{bstart + PAD}", seqb),
    ]
    # Get the aligned sequences
    bamfile = "/Users/bao/projects/female-restitution/SoSs-2024/9208.bam"
    bam = pysam.AlignmentFile(bamfile)
    chrom, start, end = a, aend - PAD, aend + PAD
    aligned_seqs = []
    for read in bam.fetch(chrom, start, end):
        try:
            read_start = read.get_reference_positions(full_length=True).index(start)
            read_end = read.get_reference_positions(full_length=True).index(end)
        except ValueError:
            continue
        aligned_sequence = read.query_sequence[read_start:read_end]
        aligned_seqs.append((read.qname, aligned_sequence))

    # Going through the reads, find the best 5 alignments to seqa and seqb
    ranked_a, ranked_b = [], []
    for read, aligned_sequence in aligned_seqs:
        alignments_a = pairwise2.align.globalxx(
            aligned_sequence, seqa, one_alignment_only=False, score_only=True
        )
        alignments_b = pairwise2.align.globalxx(
            aligned_sequence, seqb, one_alignment_only=False, score_only=True
        )
        if read in ban:
            continue
        ranked_a.append((alignments_a, read, aligned_sequence))
        ranked_b.append((alignments_b, read, aligned_sequence))
    ranked_a.sort(reverse=True)
    ranked_b.sort(reverse=True)
    SELECT = 3
    all_seqs = seqs + ranked_a[:SELECT] + ranked_b[:SELECT]
    with open(aligned_fasta, "w") as fw:
        for _, name, seq in all_seqs:
            print(">" + name, file=fw)
            print(seq, file=fw)
    logger.info("SUCCESS!")


def main(args: List[str]):
    bps = pd.read_excel("All-F1-breakpoints.xlsx")
    for _, row in bps.iterrows():
        crossover_id = row["Crossover ID"]
        if not crossover_id.startswith("9208"):
            continue
        left = row["Left"]
        right = row["Right"]
        region = f"{left} {right}"
        fasta_file = f"{crossover_id}.fasta"
        ban = set()
        if crossover_id == "9208-So-02":
            ban = set(
                [
                    "m84072_230515_103837_s3/220072159/ccs",
                    "m84072_230515_103837_s3/165938098/ccs",
                ]
            )
        write_aligned_fasta(region, ban, fasta_file)
        cmd = (
            f"clustalo -i {fasta_file} -o {fasta_file}.aln  --outfmt clu --force --full"
        )
        sh(cmd)


if __name__ == "__main__":
    main(sys.argv[1:])
