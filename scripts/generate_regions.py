import os.path as op
import sys

from collections import Counter, defaultdict
from typing import List

import pandas as pd


def load_bed(bed: str) -> pd.DataFrame:
    """
    Load a MOSDEPTH bed file into a pandas dataframe
    """
    df = pd.read_csv(bed, sep="\t")
    df.columns = ["chrom", "start", "end", "depth"]
    return df


def main(args: List[str]):
    # Find the regions that have high depth per chrom
    regions = defaultdict(list)
    if len(args) == 3:  # child, parent1, parent2
        abed, p1_bed, p2_bed = args[:3]
        af = load_bed(abed)
        p1f = load_bed(p1_bed)
        p2f = load_bed(p2_bed)

        df = af.copy()
        df["depth"] = af["depth"] / (p1f["depth"] + p2f["depth"] + 1)
    elif len(args) == 2:  # child, parents
        abed, p1_bed = args[:2]
        af = load_bed(abed)
        p1f = load_bed(p1_bed)

        df = af.copy()
        df["depth"] = af["depth"] / (p1f["depth"] + 1)

    for _, row in df.iterrows():
        chrom = row["chrom"]
        start = row["start"]
        end = row["end"]
        depth = row["depth"]
        regions[chrom].append((start, end, depth))

    d = []
    selected = []
    for chrom, data in regions.items():
        if "Chr" not in chrom:
            continue
        regions[chrom] = sorted(data, key=lambda x: x[2], reverse=True)
        chrom_selected = [
            (chrom, a, b, str(round(c)))
            for (a, b, c) in regions[chrom]
            if 5 <= c <= 100
        ]
        selected += chrom_selected
        d.append(
            (
                chrom,
                ",".join(f"{chrom}:{a}-{b}:{c}" for (chrom, a, b, c) in chrom_selected),
            )
        )

    kf = pd.DataFrame(d)
    kf.columns = ["Chrom", "Regions"]
    prefix = op.basename(abed).split(".")[0]
    poi_tsv = f"{prefix}.poi.tsv"
    kf.to_csv(poi_tsv, sep="\t", index=False)
    print(f"Points of interests written to `{poi_tsv}`")

    # Merge regions that are close to each other
    selected.sort()
    merged = []
    for i, cur in enumerate(selected):
        if i == 0:
            merged.append(cur)
            continue
        prev = merged[-1]
        if prev[0] == cur[0] and prev[2] + 20000 >= cur[1]:
            merged[-1] = (prev[0], prev[1], max(prev[2], cur[2]), f"{prev[3]},{cur[3]}")
        else:
            merged.append(cur)

    # Write the merged regions to a file
    regions_file = f"{prefix}.regions.tsv"
    counter: Counter[str] = Counter()
    with open(regions_file, "w", encoding="utf-8") as fw:
        for chrom, start, end, score in merged:
            print(f"{chrom}:{start}-{end}\t{score}", file=fw)
            counter[chrom[:2]] += 1
    print(f"Merged regions written to `{regions_file}`")
    print("Region counts:", counter)


if __name__ == "__main__":
    main(sys.argv[1:])
