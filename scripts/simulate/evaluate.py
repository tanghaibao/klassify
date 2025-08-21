import argparse
import re

import numpy as np
import pandas as pd

from pathlib import Path

from jcvi.apps.base import logger


def evaluate(data_dir: str, output: str, max_distance: int):
    """
    Combine all TSV files in the given directory into a single TSV file.
    """
    all_files = Path(data_dir).glob("*.tsv")
    results = []
    for file in all_files:
        df = pd.read_csv(file, sep="\t")
        results.append(df)
    df = pd.concat(results, ignore_index=True)
    df["Ploidy"] = df["Run"].str.count(r"\+") + 1
    df.to_csv(output, sep="\t", index=False)
    logger.info(f"Results saved to `%s`", output)

    # Process the DataFrame to find the closest true breakpoints
    matched_df = process(df)
    matched_output = output.replace(".tsv", ".matched.tsv")
    matched_df.to_csv(matched_output, sep="\t", index=False)
    logger.info(f"Matched results saved to `%s`", matched_output)

    # Print precision and recall by ploidy
    metrics = []
    for ploidy in (1, 2, 4):
        res = evaluate_recall_precision_by_ploidy(matched_df, ploidy, max_distance)
        metrics.append(res)
    metrics = pd.DataFrame(metrics)
    print(metrics)


def evaluate_recall_precision_by_ploidy(
    matched_df: pd.DataFrame, ploidy: int, max_distance: int = 50000
) -> dict:
    """
    Evaluate recall and precision by ploidy from the matched DataFrame.
    """
    kf = matched_df[matched_df["Ploidy"] == ploidy]
    recalls_df = kf[kf["Type"] == "Recall"]
    matched_recalls = recalls_df[recalls_df["Total_distance"] <= max_distance]
    precision_df = kf[kf["Type"] == "Precision"]
    matched_precision = precision_df[precision_df["Total_distance"] <= max_distance]
    recalls = len(matched_recalls) / len(recalls_df) if len(recalls_df) > 0 else 0
    precision = (
        len(matched_precision) / len(precision_df) if len(precision_df) > 0 else 0
    )
    median_distance = (
        int(matched_precision["Total_distance"].median())
        if not matched_precision.empty
        else np.nan
    )
    return {
        "Ploidy": ploidy,
        "Simulated": len(recalls_df),
        "KLASSIFY": len(precision_df),
        "Matched": len(matched_recalls),
        "Recall": recalls,
        "Precision": precision,
        "Median dist: Truth vs. KLASSIFY": median_distance,
    }


def parse_breakpoint_range(bp_str):
    """
    Parse a breakpoint string into (chrom, start, end).

    Accepted forms:
      - "chr1:12345"
      - "chr1:12345-23456"
    Returns (None, None, None) if unparsable or redacted/NaN.
    """
    if pd.isna(bp_str):
        return None, None, None
    s = str(bp_str).strip()
    if not s or "[REDACTED" in s:
        return None, None, None

    m = re.match(r"^([^:]+):\s*(\d+)(?:\s*-\s*(\d+))?$", s)
    if not m:
        return None, None, None

    chrom = m.group(1)
    start = int(m.group(2))
    end = int(m.group(3)) if m.group(3) is not None else start
    if start > end:
        start, end = end, start
    return chrom, start, end


def range_distance(s1, e1, s2, e2):
    """
    Distance between two closed intervals [s1,e1] and [s2,e2].

    - 0 if they overlap or touch
    - positive gap otherwise
    - np.inf if any endpoint is None
    """
    if any(v is None for v in (s1, e1, s2, e2)):
        return np.inf
    # overlap / touching:
    if not (e1 < s2 or e2 < s1):
        return 0
    # disjoint: gap is the min distance between edges
    if e1 < s2:
        return s2 - e1
    else:
        return s1 - e2


def find_closest_true_row(computed_row, true_rows):
    """Find the closest true row to a computed row using interval distances."""
    comp_a_chrom, comp_a_s, comp_a_e = parse_breakpoint_range(
        computed_row["A_breakpoint"]
    )
    comp_b_chrom, comp_b_s, comp_b_e = parse_breakpoint_range(
        computed_row["B_breakpoint"]
    )

    min_total_distance = np.inf
    best_match = None
    best_a_dist = None
    best_b_dist = None

    for _, true_row in true_rows.iterrows():
        true_a_chrom, true_a_s, true_a_e = parse_breakpoint_range(
            true_row["A_breakpoint"]
        )
        true_b_chrom, true_b_s, true_b_e = parse_breakpoint_range(
            true_row["B_breakpoint"]
        )

        a_dist = np.inf
        b_dist = np.inf

        if comp_a_chrom is not None and comp_a_chrom == true_a_chrom:
            a_dist = range_distance(comp_a_s, comp_a_e, true_a_s, true_a_e)

        if comp_b_chrom is not None and comp_b_chrom == true_b_chrom:
            b_dist = range_distance(comp_b_s, comp_b_e, true_b_s, true_b_e)

        # Skip if neither chromosome matches
        if a_dist == np.inf and b_dist == np.inf:
            continue

        # Sum distances for the pair (treat missing side as 0)
        total_dist = (0 if a_dist == np.inf else a_dist) + (
            0 if b_dist == np.inf else b_dist
        )

        if total_dist < min_total_distance:
            min_total_distance = total_dist
            best_match = true_row
            best_a_dist = None if a_dist == np.inf else a_dist
            best_b_dist = None if b_dist == np.inf else b_dist

    return best_match, best_a_dist, best_b_dist


def process(df: pd.DataFrame) -> pd.DataFrame:
    """
    Process the DataFrame to find the closest true breakpoints for each computed breakpoint.
    """
    # Process the data
    results = []

    for run, group in df.groupby("Run"):
        true_rows = group[group["Source"] == "true"]
        computed_rows = group[group["Source"] == "computed"]
        for _, true_row in true_rows.iterrows():
            closest_true, a_dist, b_dist = find_closest_true_row(
                true_row, computed_rows
            )
            result = {
                "Run": run,
                "Ploidy": true_row["Ploidy"],
                "Type": "Recall",
                "Computed_A_breakpoint": (
                    closest_true["A_breakpoint"] if closest_true is not None else None
                ),
                "Computed_B_breakpoint": (
                    closest_true["B_breakpoint"] if closest_true is not None else None
                ),
                "Closest_True_A_breakpoint": true_row["A_breakpoint"],
                "Closest_True_B_breakpoint": true_row["B_breakpoint"],
                "A_distance": a_dist,
                "B_distance": b_dist,
                "Total_distance": (
                    (a_dist or 0) + (b_dist or 0)
                    if (a_dist is not None or b_dist is not None)
                    else None
                ),
            }
            results.append(result)

        for _, comp_row in computed_rows.iterrows():
            closest_true, a_dist, b_dist = find_closest_true_row(comp_row, true_rows)

            result = {
                "Run": run,
                "Ploidy": comp_row["Ploidy"],
                "Type": "Precision",
                "Computed_A_breakpoint": comp_row["A_breakpoint"],
                "Computed_B_breakpoint": comp_row["B_breakpoint"],
                "Closest_True_A_breakpoint": (
                    closest_true["A_breakpoint"] if closest_true is not None else None
                ),
                "Closest_True_B_breakpoint": (
                    closest_true["B_breakpoint"] if closest_true is not None else None
                ),
                "A_distance": a_dist,
                "B_distance": b_dist,
                "Total_distance": (
                    (a_dist or 0) + (b_dist or 0)
                    if (a_dist is not None or b_dist is not None)
                    else None
                ),
            }
            results.append(result)

    result_df = pd.DataFrame(results)
    result_df["A_distance"] = result_df["A_distance"].astype("Int64")
    result_df["B_distance"] = result_df["B_distance"].astype("Int64")
    result_df["Total_distance"] = result_df["Total_distance"].astype("Int64")
    return result_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate a model on a dataset.")
    parser.add_argument("data_dir", type=str, help="Path to the model to evaluate.")
    parser.add_argument(
        "--output",
        type=str,
        default="results.tsv.gz",
        help="File to save evaluation results.",
    )
    parser.add_argument(
        "--max-distance",
        type=int,
        default=2000,
        help="Max distance for matching breakpoints.",
    )

    args = parser.parse_args()

    evaluate(args.data_dir, args.output, args.max_distance)
