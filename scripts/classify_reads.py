import sys

from glob import glob
from multiprocessing import Pool
from typing import List, Set, Tuple

import pandas as pd


def get_pairs() -> Set[Tuple[str, str]]:
    df = pd.read_csv("chimeric_chrom_pairs.tsv", sep="\t")
    pairs = set()
    for _, row in df.iterrows():
        a = row["ChromA"]
        b = row["ChromB"]
        pairs.add(tuple(sorted((a, b))))
    return pairs


def get_reads(rc: str) -> pd.DataFrame:
    rf = pd.read_csv(rc, sep="\t")
    filtered = []
    for _, row in rf.iterrows():
        kmers = row["Kmers"]
        classification = row["Classification"]
        if "Unclassified" in classification:
            continue
        ab, scores = classification.split(":", 1)
        a, b = ab.replace(".fa", "").split(",")
        if a[:7] != b[:7]:
            continue
        a, b = tuple(sorted((a, b)))
        ascore, bscore = scores.split(",", 1)
        ascore, bscore = int(ascore), int(bscore)
        if kmers >= 300 and ascore + bscore >= 50 and bscore >= 10:
            row["Label"] = f"{a}_{b}"
            filtered.append(row)
    return pd.DataFrame(filtered)


def main(args: List[str]):
    workdir = args[0]
    rcs = glob(f"{workdir}/*.read_classifications.tsv")
    print(f"Found f{len(rcs)} read classifications")
    pool = Pool(100)
    dfs = pool.map(get_reads, rcs)
    df = pd.concat(dfs)
    df.to_csv(f"{workdir}.filtered.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main(sys.argv[1:])
