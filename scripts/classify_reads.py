from glob import glob
from multiprocessing import Pool
from typing import Set, Tuple

import pandas as pd


def get_pairs() -> Set[Tuple[str, str]]:
    df = pd.read_csv("chimeric_chrom_pairs.tsv", sep="\t")
    pairs = set()
    for _, row in df.iterrows():
        a = row["ChromA"]
        b = row["ChromB"]
        pairs.add(tuple(sorted((a, b))))
    return pairs


def get_reads(args) -> pd.DataFrame:
    rc, pairs = args
    rf = pd.read_csv(rc, sep="\t")
    filtered = []
    for _, row in rf.iterrows():
        length = row["Length"]
        kmers = row["Kmers"]
        classification = row["Classification"]
        if "Unclassified" in classification:
            continue
        ab, scores = classification.split(":", 1)
        a, b = ab.replace(".fa", "").split(",")
        a, b = tuple(sorted((a, b)))
        ascore, bscore = scores.split(",", 1)
        ascore, bscore = int(ascore), int(bscore)
        if kmers >= 100 and length > 1000 and ascore + bscore >= 50 and bscore >= 10 and (a, b) in pairs:
            row["Label"] = f"{a}_{b}"
            filtered.append(row)
    return pd.DataFrame(filtered)


def main():
    pairs = get_pairs()
    print(pairs)
    rcs = glob("LAp_hifi_classify/*.read_classifications.tsv")
    args = [(rc, pairs) for rc in rcs]
    pool = Pool(100)
    dfs = pool.map(get_reads, args)
    df = pd.concat(dfs)
    df.to_csv("filtered_LAP_HIFI.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
