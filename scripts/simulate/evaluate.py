import argparse

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
    return {
        "Ploidy": ploidy,
        "Simulated": len(recalls_df),
        "KLASSIFY": len(precision_df),
        "Matched": len(matched_recalls),
        "Recall": recalls,
        "Precision": precision,
    }


def parse_breakpoint(bp_str):
    """Parse breakpoint string to extract chromosome and position"""
    if pd.isna(bp_str) or "[REDACTED" in str(bp_str):
        return None, None

    try:
        parts = str(bp_str).split(":")
        if len(parts) != 2:
            return None, None

        chrom = parts[0]
        pos_part = parts[1]

        # Handle range format (start-end)
        if "-" in pos_part:
            start, end = map(int, pos_part.split("-"))
            midpoint = (start + end) / 2
            return chrom, midpoint
        else:
            pos = int(pos_part)
            return chrom, pos
    except:
        return None, None


def calculate_distance(pos1, pos2):
    """Calculate distance between two positions"""
    if pos1 is None or pos2 is None:
        return np.inf
    return abs(pos1 - pos2)


def find_closest_true_row(computed_row, true_rows):
    """Find the closest true row to a computed row"""
    comp_a_chrom, comp_a_pos = parse_breakpoint(computed_row["A_breakpoint"])
    comp_b_chrom, comp_b_pos = parse_breakpoint(computed_row["B_breakpoint"])

    min_total_distance = np.inf
    best_match = None
    best_a_dist = None
    best_b_dist = None

    for _, true_row in true_rows.iterrows():
        true_a_chrom, true_a_pos = parse_breakpoint(true_row["A_breakpoint"])
        true_b_chrom, true_b_pos = parse_breakpoint(true_row["B_breakpoint"])

        # Calculate distances only for matching chromosomes
        a_dist = np.inf
        b_dist = np.inf

        if comp_a_chrom == true_a_chrom:
            a_dist = calculate_distance(comp_a_pos, true_a_pos)

        if comp_b_chrom == true_b_chrom:
            b_dist = calculate_distance(comp_b_pos, true_b_pos)

        # Skip if no matching chromosomes
        if a_dist == np.inf and b_dist == np.inf:
            continue

        # Calculate total distance (you can modify this logic as needed)
        # Here I'm using the sum of distances, but you could use max, euclidean, etc.
        total_dist = (a_dist if a_dist != np.inf else 0) + (
            b_dist if b_dist != np.inf else 0
        )

        if total_dist < min_total_distance:
            min_total_distance = total_dist
            best_match = true_row
            best_a_dist = a_dist if a_dist != np.inf else None
            best_b_dist = b_dist if b_dist != np.inf else None

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
        default=20000,
        help="Max distance for matching breakpoints.",
    )

    args = parser.parse_args()

    evaluate(args.data_dir, args.output, args.max_distance)
