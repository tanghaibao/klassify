import sys

from dataclasses import dataclass
from typing import List

from jcvi.apps.base import logger, sh


FLANKING = 10000


@dataclass
class Region:
    """
    A region in the genome.
    """

    chrom: str
    start: int
    end: int
    name: str

    @classmethod
    def from_str(cls, region: str) -> "Region":
        """
        Create a Region object from a string.
        """
        chrom, start_end = region.split(":")
        start, end = start_end.split("-")
        return cls(chrom, int(start), int(end), region)

    def __str__(self):
        return f"{self.chrom}:{self.start}-{self.end}"


def main(regions_file: str, genome_fasta: str, *bam_files: List[str]):
    """
    Generate an IGV report from a list of regions.
    """
    regions = []
    with open(regions_file, "r", encoding="utf-8") as f:
        for row in f:
            region, _ = row.strip().split()
            regions.append(Region.from_str(region))

    # Make a BED file
    bed_file = regions_file.rsplit(".", 1)[0] + ".bed"
    with open(bed_file, "w", encoding="utf-8") as f:
        for region in regions:
            f.write(f"{region.chrom}\t{region.start}\t{region.end}\t{region.name}\n")
    logger.info("Generated BED file: %s", bed_file)

    # Call create_report
    logger.info("Creating report for %s...", regions_file)
    html_file = regions_file.rsplit(".", 1)[0] + ".html"
    cmd = f"create_report {bed_file}"
    cmd += f" --fasta {genome_fasta}"
    cmd += f" --flanking {FLANKING}"
    cmd += f" --tracks {' '.join(bam_files)}"
    cmd += f" --output {html_file}"
    sh(cmd)

    logger.info("Generated HTML report: %s", html_file)
    logger.info("Success")


if __name__ == "__main__":
    args = sys.argv[1:]
    main(*args)
