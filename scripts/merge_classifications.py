import sys

from glob import glob
from typing import List

import pandas as pd


def main(args: List[str]):
    workdir = args[0]
    flist = glob(f"{workdir}/*.read_classifications.tsv")
    data = []
    for f in flist:
        df = pd.read_csv(f, sep="\t")
        data.append(df)
    merged = pd.concat(data)
    merged.sort_values("ID", inplace=True)
    out_tsv = f"{workdir}.filtered.tsv"
    merged.to_csv(out_tsv, sep="\t", index=False)
    print(f"Classification merged and written to `{out_tsv}`")


if __name__ == "__main__":
    main(sys.argv[1:])
