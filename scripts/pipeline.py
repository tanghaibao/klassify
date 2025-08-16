#!/usr/bin/env python3
"""
Orchestrate the KLASSIFY-based pipeline:

Inputs:
  - f1_reads.fa
  - parent_reads.fa
  - parents.genome.fa

External tools expected in PATH:
  - klassify
  - minimap2
  - samtools
  - bamToBed   (bedtools bamtobed; alias often installed as `bamToBed`)

Two helper scripts (provide via --scripts-dir or absolute paths):
  - split_reads_by_breakpoint.py
  - cluster_to_paired_regions.py
"""

from __future__ import annotations
import argparse
from pathlib import Path
import sys

from simulate.utils import check_tool, run, run_pipe
from jcvi.apps.base import logger, mkdir, need_update

DEFAULT_SCRIPT_PATH = Path(__file__).parent


def main():
    ap = argparse.ArgumentParser(description="KLASSIFY pipeline (Python wrapper).")
    ap.add_argument("f1_reads", type=Path, help="F1 reads FASTA")
    ap.add_argument("parent_reads", type=Path, help="Parental reads FASTA")
    ap.add_argument("parents_ref", type=Path, help="Combined parental reference FASTA")
    ap.add_argument(
        "--workdir", type=Path, default=Path("."), help="Working directory (default: .)"
    )
    ap.add_argument(
        "--scripts-dir",
        type=Path,
        default=DEFAULT_SCRIPT_PATH,
        help="Directory containing helper scripts (split_reads_by_breakpoint.py, cluster_to_paired_regions.py).",
    )
    ap.add_argument(
        "--kmers",
        type=Path,
        default=None,
        help="Prebuilt kmers.bc (skip build if provided).",
    )
    ap.add_argument(
        "--threads-minimap2", type=int, default=32, help="Threads for minimap2 (-t)"
    )
    ap.add_argument(
        "--threads-sort", type=int, default=8, help="Threads for samtools sort (-@)"
    )
    ap.add_argument(
        "--preset",
        default="map-hifi",
        choices=["map-hifi", "map-ont"],
        help="Minimap2 preset (default: map-hifi)",
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="Overwrite/redo steps even if outputs exist",
    )
    ap.add_argument(
        "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"]
    )
    args = ap.parse_args()

    logger.setLevel(args.log_level)

    input_files = (args.f1_reads, args.parent_reads, args.parents_ref)

    # sanity checks
    for p in input_files:
        if not p.exists():
            ap.error(f"Input not found: {p}")

    for tool in ("klassify", "minimap2", "samtools", "bamToBed"):
        check_tool(tool)

    work = args.workdir.resolve()
    mkdir(work)

    # paths
    f1_out_dir = work / "f1_classify"
    parent_out_dir = work / "parent_classify"
    mkdir(f1_out_dir)
    mkdir(parent_out_dir)

    kmers_bc = args.kmers.resolve() if args.kmers else (work / "kmers.bc")

    f1_filtered_tsv = work / "f1_classify.filtered.tsv"
    parent_filtered_tsv = work / "parent_classify.filtered.tsv"

    f1_extracted_fa = work / "f1_classify.fa"
    parent_extracted_fa = work / "parent_classify.fa"

    f1_bam = work / "f1_classify.bam"
    parent_bam = work / "parent_classify.bam"

    regions_tsv = work / "f1_classify.regions.tsv"
    regions_fa = work / "f1_classify.regions.fasta"
    roi_bam = work / "f1_classify.roi.bam"
    roi_bed = work / "f1_classify.roi.bed"
    roi_tsv = work / "f1_classify.roi.tsv"
    paired_regions = work / "f1_classify.paired.regions"

    # helper scripts
    scripts_dir = args.scripts_dir.resolve() if args.scripts_dir else None
    split_script = (
        (scripts_dir / "split_reads_by_breakpoint.py")
        if scripts_dir
        else Path("split_reads_by_breakpoint.py")
    )
    cluster_script = (
        (scripts_dir / "cluster_to_paired_regions.py")
        if scripts_dir
        else Path("cluster_to_paired_regions.py")
    )

    if not split_script.exists():
        logger.warning(
            "split_reads_by_breakpoint.py not found at %s (will try anyway if it's in PATH).",
            split_script,
        )
    if not cluster_script.exists():
        logger.warning(
            "cluster_to_paired_regions.py not found at %s (will try anyway if it's in PATH).",
            cluster_script,
        )

    # 1) Build unique k-mers from parental genomes
    if args.kmers and kmers_bc.exists():
        logger.info("Using prebuilt k-mer DB: %s", kmers_bc)
    elif not need_update(args.parents_ref, kmers_bc) and not args.force:
        logger.info("Found kmers DB, skipping build: %s", kmers_bc)
    else:
        run(["klassify", "build", str(args.parents_ref), "-o", str(kmers_bc)], cwd=work)

    # 2) Classify F1, extract chimeric reads, map to parents reference
    if not need_update([kmers_bc, args.f1_reads], f1_filtered_tsv) and not args.force:
        logger.info("F1 classify already done: %s", f1_filtered_tsv)
    else:
        run(
            [
                "klassify",
                "classify",
                str(kmers_bc),
                str(args.f1_reads),
                "-o",
                "f1_classify",
            ],
            cwd=work,
        )

    if (
        not need_update([f1_filtered_tsv, args.f1_reads], f1_extracted_fa)
        and not args.force
    ):
        logger.info("F1 extract already done: %s", f1_extracted_fa)
    else:
        run(
            [
                "klassify",
                "extract",
                str(f1_filtered_tsv),
                str(args.f1_reads),
                "-o",
                str(f1_extracted_fa),
            ],
            cwd=work,
        )

    if not need_update([f1_extracted_fa, args.parents_ref], f1_bam) and not args.force:
        logger.info("F1 alignment already exists: %s", f1_bam)
    else:
        run_pipe(
            [
                "minimap2",
                "-t",
                str(args.threads_minimap2),
                "-ax",
                args.preset,
                "--eqx",
                "--secondary=no",
                str(args.parents_ref),
                str(f1_extracted_fa),
            ],
            ["samtools", "sort", "-@", str(args.threads_sort), "-o", str(f1_bam)],
            cwd=work,
            stdout_path=f1_bam,
        )

    # 3) Classify parent reads (control), extract, map
    if (
        not need_update([kmers_bc, args.parent_reads], parent_filtered_tsv)
        and not args.force
    ):
        logger.info("Parent classify already done: %s", parent_filtered_tsv)
    else:
        run(
            [
                "klassify",
                "classify",
                str(kmers_bc),
                str(args.parent_reads),
                "-o",
                "parent_classify",
            ],
            cwd=work,
        )

    if (
        not need_update([parent_filtered_tsv, args.parent_reads], parent_extracted_fa)
        and not args.force
    ):
        logger.info("Parent extract already done: %s", parent_extracted_fa)
    else:
        run(
            [
                "klassify",
                "extract",
                str(parent_filtered_tsv),
                str(args.parent_reads),
                "-o",
                str(parent_extracted_fa),
            ],
            cwd=work,
        )

    if (
        not need_update([parent_extracted_fa, args.parents_ref], parent_bam)
        and not args.force
    ):
        logger.info("Parent alignment already exists: %s", parent_bam)
    else:
        run_pipe(
            [
                "minimap2",
                "-t",
                str(args.threads_minimap2),
                "-ax",
                args.preset,
                "--eqx",
                "--secondary=no",
                str(args.parents_ref),
                str(parent_extracted_fa),
            ],
            ["samtools", "sort", "-@", str(args.threads_sort), "-o", str(parent_bam)],
            cwd=work,
            stdout_path=parent_bam,
        )

    # 4) Regions present in F1 but NOT parent (control)
    if not need_update([f1_bam, parent_bam], regions_tsv) and not args.force:
        logger.info("Regions TSV already exists: %s", regions_tsv)
    else:
        run(["klassify", "regions", str(f1_bam), str(parent_bam)], cwd=work)

    # 5) Refine the breakpoints
    # Extract BAM segments overlapping regions
    if not need_update([regions_tsv, f1_bam], regions_fa) and not args.force:
        logger.info("Regions FASTA already exists: %s", regions_fa)
    else:
        run(["klassify", "extract-bam", str(regions_tsv), str(f1_bam)], cwd=work)
        # Breakpoint calling on extracted reads
        # (output files are produced next to input; klassify decides names internally)
        run(["klassify", "breakpoint", str(kmers_bc), str(regions_fa)], cwd=work)

    # Split reads at breakpoints
    split_input = regions_fa
    split_output = work / "f1_classify.regions.split.fasta"
    if not need_update(split_input, split_output) and not args.force:
        logger.info("Split reads already exist: %s", split_output)
    else:
        run([sys.executable, str(split_script), str(split_input)], cwd=work)

    # Remap the split reads to parents to get a crisp ROI BAM
    if not need_update(split_output, roi_bam) and not args.force:
        logger.info("ROI BAM already exists: %s", roi_bam)
    else:
        run_pipe(
            [
                "minimap2",
                "-t",
                str(args.threads_minimap2),
                "-ax",
                args.preset,
                "--eqx",
                "--secondary=no",
                str(args.parents_ref),
                str(split_output),
            ],
            ["samtools", "sort", "-@", str(args.threads_sort), "-o", str(roi_bam)],
            cwd=work,
            stdout_path=roi_bam,
        )

    # BAM → BED
    if not need_update(roi_bam, roi_bed) and not args.force:
        logger.info("ROI BED already exists: %s", roi_bed)
    else:
        with open(roi_bed, "wb") as out:
            run(["bamToBed", "-i", str(roi_bam)], cwd=work, stdout=out)

    # Cluster paired regions
    if not need_update(roi_bed, paired_regions) and not args.force:
        logger.info("Paired regions TSV already exists: %s", paired_regions)
    else:
        with open(roi_tsv, "wb") as out:
            run(
                [sys.executable, str(cluster_script), str(roi_bed)],
                cwd=work,
                stdout=out,
            )

    logger.info("Done.")
    logger.info(
        "Key outputs:\n"
        "  k-mers:          %s\n"
        "  F1 BAM:          %s\n"
        "  Parent BAM:      %s\n"
        "  Regions TSV:     %s\n"
        "  Regions FASTA:   %s\n"
        "  ROI BAM:         %s\n"
        "  ROI BED:         %s\n"
        "  ROI TSV:         %s\n"
        "  Paired Regions:  %s",
        kmers_bc,
        f1_bam,
        parent_bam,
        regions_tsv,
        regions_fa,
        roi_bam,
        roi_bed,
        roi_tsv,
        paired_regions,
    )


if __name__ == "__main__":
    main()
