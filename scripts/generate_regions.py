from collections import defaultdict
from typing import List

import pandas as pd
import sys


def load_bed(bed: str) -> pd.DataFrame:
    """
    Load a MOSDEPTH bed file into a pandas dataframe
    """
    df = pd.read_csv(bed, sep="\t")
    df.columns = ["chrom", "start", "end", "depth"]
    return df


def main(args: List[str]):
    abed, bed = args[:2]
    hf = load_bed(abed)
    of = load_bed(bed)

    # Find the regions that have high depth per chrom
    regions = defaultdict(list)
    df = hf.copy()
    df["depth"] = hf["depth"] / (of["depth"] + 1)
    for _, row in df.iterrows():
        chrom = row["chrom"]
        start = row["start"]
        end = row["end"]
        depth = row["depth"]
        regions[chrom].append((start, end, depth))

    d = []
    for chrom, data in regions.items():
        if "utg" in chrom:
            continue
        regions[chrom] = sorted(data, key=lambda x: x[2], reverse=True)
        d.append(
            (
                chrom,
                ",".join(
                    f"{chrom}:{a}-{b}:{c}" for (a, b, c) in regions[chrom] if c > 5
                ),
            )
        )

    kf = pd.DataFrame(d)
    kf.columns = ["Chrom", "Regions"]
    kf.to_csv("poi.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main(sys.argv[1:])
