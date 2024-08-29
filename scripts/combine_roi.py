import os.path as op
import sys

from typing import List

import pandas as pd

from jcvi.apps.base import logger


def main(args: List[str]):
    if len(args) < 2:
        print("Usage: python combine_roi.py <roi1.bed> <roi2.bed> ...")
        sys.exit(1)

    logger.info(f"Combining %d ROIs", len(args))

    rois = []
    for arg in args:
        if arg.startswith("combined"):
            logger.info("Skipping `%s`", arg)
            continue
        df = pd.read_csv(arg, sep="\t")
        pid = op.basename(arg).split("-")[0]
        df["Crossover ID"] = df["Crossover ID"].apply(lambda x: f"{pid}-{x}")
        rois.append(df)

    df = pd.concat(rois)
    df.to_csv("combined-roi-with-reads.tsv", sep="\t", header=True, index=False)
    df_no_reads = df.dropna(subset=["Left"])
    combined_roi = "combined-roi.tsv"
    df_no_reads.to_csv(combined_roi, sep="\t", header=True, index=False)
    logger.info("Combined ROIs written to `%s`", combined_roi)


if __name__ == "__main__":
    main(sys.argv[1:])
