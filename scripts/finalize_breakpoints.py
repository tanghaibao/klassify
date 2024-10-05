import sys

from collections import Counter
from typing import List

import pandas as pd

from jcvi.apps.base import logger


def main(args: List[str]):
    if len(args) != 1:
        print("Usage: finalize_breakpoints.py All-F1-breakpoints-with-reads-IGV.xlsx")
        sys.exit(1)

    (proofed_xls,) = args
    df = pd.read_excel(proofed_xls, sheet_name="combined-roi-with-reads")
    counter = Counter()
    summary = Counter()
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
        accession, _ = crossover_id.split("-", 1)
        if accession.startswith("92") and gamete == "So" and left_type == right_type:
            removed["So-Type match"] += 1
            continue
        if gamete == "Ss" and "Type II" in (left_type, right_type):
            removed["Ss-Type mismatch"] += 1
            continue
        counter[(accession, gamete)] += 1
        # summary[(accession, gamete, "+".join(sorted([left_type,
        # right_type])))] += 1
        summary[(accession, gamete, left_type)] += 1
        summary[(accession, gamete, right_type)] += 1
        new_crossover_id = f"{accession}-{gamete}-{counter[(accession, gamete)]:02d}"
        old_to_new[crossover_id] = new_crossover_id
    logger.info("Removed: %s", removed)

    # Update the Crossover ID
    df["Crossover ID"] = df["Crossover ID"].map(old_to_new)
    df = df[~df["Crossover ID"].isna()]
    new_xls = proofed_xls.replace("-IGV.xlsx", ".xlsx")
    df.to_excel(new_xls, index=False)
    logger.info("Updated file saved to `%s`", new_xls)

    # Save a summary table
    summary_df = []
    accessions = sorted(set(accession for accession, _ in counter))
    for accession in accessions:
        summary_df.append(
            {
                "F1": accession,
                # "So Type I+Type II": summary[(accession, "So", "Type I+Type II")],
                # "Ss Type I+Type I": summary[(accession, "Ss", "Type I+Type
                # I")],
                "So Type I": summary[(accession, "So", "Type I")],
                "So Type II": summary[(accession, "So", "Type II")],
                "Ss Type I": summary[(accession, "Ss", "Type I")],
                "Ss Type II": summary[(accession, "Ss", "Type II")],
                "Total Pairs": (
                    counter[(accession, "So")] + counter[(accession, "Ss")]
                ),
            }
        )
    summary_df = pd.DataFrame(summary_df)
    print(summary_df)
    summary_xls = "F1-breakpoints-summary.xlsx"
    summary_df.to_excel(summary_xls, index=False)
    logger.info("Summary table saved to `%s`", summary_xls)

    # Write a simplified crossover table without the reads
    df = df[["Crossover ID", "Left", "Right", "Left Type", "Right Type"]]
    df.dropna(subset=["Left", "Right"], inplace=True)
    simple_xls = new_xls.replace("-with-reads", "")
    df.to_excel(simple_xls, index=False)
    logger.info("Simplifed crossover table saved to `%s`", simple_xls)

    # make_roi.py
    roi = []
    for _, row in df.iterrows():
        accession, _ = row["Crossover ID"].split("-", 1)
        accession_bed_gz = f"{accession}Aln2ParentsHiFi_2024.regions.bed.gz"
        roi += [
            {
                "Accession": accession_bed_gz,
                "Region": row["Left"],
                "Type": row["Left Type"].replace("Type ", ""),
            },
            {
                "Accession": accession_bed_gz,
                "Region": row["Right"],
                "Type": row["Right Type"].replace("Type ", ""),
            },
        ]
    roi_df = pd.DataFrame(roi)
    roi_csv = "roi.csv"
    roi_df.to_csv(roi_csv, index=False, header=None)
    logger.info("ROI table saved to `%s`", roi_csv)


if __name__ == "__main__":
    main(sys.argv[1:])
