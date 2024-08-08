import json
import os.path as op
import sys

from dataclasses import dataclass
from typing import List

from jcvi.apps.base import logger, sh


FLANKING = 50000
TRACK_CONFIG = {
    "type": "alignment",
    "hideSmallIndels": True,
    "indelSizeThreshold": 2,
    "height": 250,
    # "displayMode": "SQUISHED",
}


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
        return cls(chrom, int(start) - 1, int(end), region)


def main(regions_file: str, genome_fasta: str, *bam_files: List[str]):
    """
    Generate an IGV report from a list of regions.
    """
    prefix = regions_file.rsplit(".", 1)[0]
    regions = []
    with open(regions_file, "r", encoding="utf-8") as f:
        for row in f:
            region, _ = row.strip().split()
            regions.append(Region.from_str(region))

    # Make a BED file
    bed_file = prefix + ".bed"
    with open(bed_file, "w", encoding="utf-8") as f:
        for region in regions:
            f.write(f"{region.chrom}\t{region.start}\t{region.end}\t{region.name}\n")
    logger.info("Generated BED file: %s", bed_file)

    # Write track configuration
    track_config = prefix + ".track_config.json"
    all_tracks = []
    for bam_file in bam_files:
        track = TRACK_CONFIG.copy()
        track["name"] = op.basename(bam_file)
        track["url"] = bam_file
        all_tracks.append(track)
    with open(track_config, "w", encoding="utf-8") as f:
        f.write(json.dumps(all_tracks, indent=4))
    logger.info("Generated track configuration: %s", track_config)

    # Call create_report
    logger.info("Creating report for %s...", regions_file)
    html_file = prefix + ".html"
    cmd = f"create_report {bed_file}"
    cmd += f" --fasta {genome_fasta}"
    cmd += f" --flanking {FLANKING}"
    cmd += f" --track-config {track_config}"
    cmd += f" --output {html_file}"
    sh(cmd)

    logger.info("Generated HTML report: %s", html_file)
    logger.info("Success")


if __name__ == "__main__":
    args = sys.argv[1:]
    main(*args)
