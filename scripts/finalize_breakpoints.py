import sys

from collections import Counter
from typing import List

import pandas as pd

from jcvi.apps.base import logger


def main(args: List[str]):
    if len(args) != 1:
        print("Usage: finalize_breakpoints.py All-F1-Breakpoints-with-reads.xlsx")
        sys.exit(1)

    (proofed_xls,) = args
    df = pd.read_excel(proofed_xls, sheet_name="combined-roi-with-reads")
    counter = Counter()
    old_to_new = {}
    removed = Counter()
    for _, row in df.iterrows():
        crossover_id = row["Crossover ID"]
        left_type = row["Left Type"]
        right_type = row["Right Type"]
        if "bad" in (left_type, right_type):
            removed["bad"] += 1
            continue
        left = row["Left"]
        if pd.isna(left):
            continue
        gamete = left[:2]
        if gamete == "Ss" and "Type II" in (left_type, right_type):
            removed["Ss-Type II"] += 1
            continue
        accession, _ = crossover_id.split("-", 1)
        counter[(accession, gamete)] += 1
        new_crossover_id = f"{accession}-{gamete}-{counter[(accession, gamete)]:02d}"
        old_to_new[crossover_id] = new_crossover_id
    logger.info("Removed: %s", removed)

    # Update the Crossover ID
    df["Crossover ID"] = df["Crossover ID"].map(old_to_new)
    df = df[~df["Crossover ID"].isna()]
    new_xls = proofed_xls.replace(".xlsx", "-final.xlsx")
    df.to_excel(new_xls, index=False)
    logger.info("Updated file saved to %s", new_xls)

    # Save a summary table


if __name__ == "__main__":
    main(sys.argv[1:])
