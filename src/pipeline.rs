//! Orchestrate the end-to-end breakpoint→crossover pipeline.
//!
//! Purpose: Drive all stages from raw reads and parental references to a final,
//! high-confidence set of crossover pairs, wiring together internal modules and
//! external tools with sensible defaults, parallelism, and resumability.
//!
//! Typical inputs:
//! - F1 reads (FA/FASTQ)
//! - Parental read sets and/or combined parental reference (FASTA)
//! - Parameters: k, window size, clustering radius, distance limits, min supports,
//!   thread count, temp/output directories, and flags like --resume/--keep-temp.
//!
//! Main stages (enabled by flags and auto-skipped if outputs exist):
//! 1) Preflight: validate inputs, check external tools in PATH, create work dirs.
//! 2) Index: build or load parent-specific unique k-mer index.
//! 3) Label: annotate F1 reads with parental origin per k-mer (streaming).
//! 4) Breakpoint calling: invoke `breakpoint` module to emit per-read calls/TSV.
//! 5) (Optional) Split: cut reads at calls to segments for mapping/debugging.
//! 6) Mapping: align segments or full reads to parental refs (e.g. minimap2),
//!    produce PAF/BAM for coordinate normalization and QC.
//! 7) Site clustering: group breakpoint calls into left/right sites with support.
//! 8) Pairing: form candidate left↔right pairs and select non-overlapping pairs
//!    with a greedy chooser; write final crossover table.
//! 9) Reports: dump TSV/BED, basic metrics, and logs for each stage.
//!
//! Outputs:
//! - Breakpoint TSV per read, site BED/TSV, paired-region (crossover) table,
//!   optional segment FASTA/FASTQ and alignment files for inspection.
//!
//! Execution model:
//! - Parallel where safe (I/O bounded stages are batched; CPU heavy stages use
//!   `--jobs`). Each stage declares its products and can be resumed via stamps.
//! - Errors fail fast with context; partial results are left in temp for audit.
//!
//! Notes:
//! - Coordinates can be normalized to a chosen reference if alignments are present;
//!   otherwise the pipeline operates in read space.
//! - Parameters balance recall vs precision; presets are provided for long reads.
//!
//! Complexity: dominated by mapping and per-read scanning; other stages are near-linear
//! in the number of breakpoints and candidate pairs.

use anyhow::{anyhow, bail, Context, Result};
use clap::{Parser, ValueEnum};
use log::{debug, info};
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::SystemTime;
use which::which;

/// KLASSIFY pipeline (Rust wrapper) — ported from pipeline.py
#[derive(Parser, Debug)]
pub struct PipelineArgs {
    /// F1 reads FASTA
    f1_reads: PathBuf,
    /// Parental reads FASTA
    parent_reads: PathBuf,
    /// Combined parental reference FASTA
    parents_ref: PathBuf,

    /// Working directory (default: .)
    #[arg(long, default_value = ".")]
    workdir: PathBuf,

    /// Threads for minimap2 (-t)
    #[arg(long, default_value_t = 32)]
    threads_minimap2: usize,

    /// Threads for samtools sort (-@)
    #[arg(long, default_value_t = 8)]
    threads_sort: usize,

    /// Minimap2 preset
    #[arg(long, value_enum, default_value_t = Preset::MapHifi)]
    preset: Preset,

    /// Overwrite/redo steps even if outputs exist
    #[arg(long)]
    force: bool,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum Preset {
    #[value(name = "map-hifi")]
    MapHifi,
    #[value(name = "map-ont")]
    MapOnt,
}

impl Preset {
    fn as_str(&self) -> &'static str {
        match self {
            Preset::MapHifi => "map-hifi",
            Preset::MapOnt => "map-ont",
        }
    }
}

/// Run the KLASSIFY pipeline. This is useful for small runs. For larger datasets, break the steps
/// into separate commands and use faSplit to split the reads into smaller chunks to increase level
/// of parallelism.
pub fn pipeline(args: PipelineArgs) -> Result<()> {
    // sanity checks
    for p in [&args.f1_reads, &args.parent_reads, &args.parents_ref] {
        if !p.exists() {
            bail!("Input not found: {}", p.display());
        }
    }

    for tool in ["klassify", "minimap2", "samtools"] {
        check_tool(tool)?;
    }

    let work = fs::canonicalize(&args.workdir).unwrap_or_else(|_| args.workdir.clone());
    mkdir(&work)?;

    // paths
    let f1_out_dir = work.join("f1_classify");
    let parent_out_dir = work.join("parent_classify");
    mkdir(&f1_out_dir)?;
    mkdir(&parent_out_dir)?;

    let kmers_bc = work.join("kmers.bc");

    let f1_filtered_tsv = work.join("f1_classify.filtered.tsv");
    let parent_filtered_tsv = work.join("parent_classify.filtered.tsv");

    let f1_extracted_fa = work.join("f1_classify.fa");
    let parent_extracted_fa = work.join("parent_classify.fa");

    let f1_bam = work.join("f1_classify.bam");
    let parent_bam = work.join("parent_classify.bam");

    let regions_tsv = work.join("f1_classify.regions.tsv");
    let regions_fa = work.join("f1_classify.regions.fasta");
    let roi_bam = work.join("f1_classify.roi.bam");
    let roi_tsv = work.join("f1_classify.roi.tsv");
    let paired_regions = work.join("f1_classify.roi.paired.regions"); // kept to mirror Python

    // 1) Build unique k-mers from parental genomes
    if up_to_date(&[&args.parents_ref], &[&kmers_bc]) && !args.force {
        info!("Found kmers DB, skipping build: {}", kmers_bc.display());
    } else {
        run(
            &work,
            "klassify",
            &[
                OsStr::new("build"),
                args.parents_ref.as_os_str(),
                OsStr::new("-o"),
                kmers_bc.as_os_str(),
            ],
        )?;
    }

    // 2) Classify F1, extract chimeric reads, map to parents reference
    if up_to_date(&[&kmers_bc, &args.f1_reads], &[&f1_filtered_tsv]) && !args.force {
        info!("F1 classify already done: {}", f1_filtered_tsv.display());
    } else {
        run(
            &work,
            "klassify",
            &[
                OsStr::new("classify"),
                kmers_bc.as_os_str(),
                args.f1_reads.as_os_str(),
                OsStr::new("-o"),
                OsStr::new("f1_classify"),
            ],
        )?;
    }

    if up_to_date(&[&f1_filtered_tsv, &args.f1_reads], &[&f1_extracted_fa]) && !args.force {
        info!("F1 extract already done: {}", f1_extracted_fa.display());
    } else {
        run(
            &work,
            "klassify",
            &[
                OsStr::new("extract"),
                f1_filtered_tsv.as_os_str(),
                args.f1_reads.as_os_str(),
                OsStr::new("-o"),
                f1_extracted_fa.as_os_str(),
            ],
        )?;
    }

    if up_to_date(&[&f1_extracted_fa, &args.parents_ref], &[&f1_bam]) && !args.force {
        info!("F1 alignment already exists: {}", f1_bam.display());
    } else {
        run_pipe(
            &work,
            (
                "minimap2",
                vec![
                    OsStr::new("-t"),
                    OsStr::new(&args.threads_minimap2.to_string()),
                    OsStr::new("-ax"),
                    OsStr::new(args.preset.as_str()),
                    OsStr::new("--eqx"),
                    OsStr::new("--secondary=no"),
                    args.parents_ref.as_os_str(),
                    f1_extracted_fa.as_os_str(),
                ],
            ),
            (
                "samtools",
                vec![
                    OsStr::new("sort"),
                    OsStr::new("-@"),
                    OsStr::new(&args.threads_sort.to_string()),
                    OsStr::new("-o"),
                    f1_bam.as_os_str(),
                ],
            ),
        )?;
    }

    // 3) Classify parent reads (control), extract, map
    if up_to_date(&[&kmers_bc, &args.parent_reads], &[&parent_filtered_tsv]) && !args.force {
        info!(
            "Parent classify already done: {}",
            parent_filtered_tsv.display()
        );
    } else {
        run(
            &work,
            "klassify",
            &[
                OsStr::new("classify"),
                kmers_bc.as_os_str(),
                args.parent_reads.as_os_str(),
                OsStr::new("-o"),
                OsStr::new("parent_classify"),
            ],
        )?;
    }

    if up_to_date(
        &[&parent_filtered_tsv, &args.parent_reads],
        &[&parent_extracted_fa],
    ) && !args.force
    {
        info!(
            "Parent extract already done: {}",
            parent_extracted_fa.display()
        );
    } else {
        run(
            &work,
            "klassify",
            &[
                OsStr::new("extract"),
                parent_filtered_tsv.as_os_str(),
                args.parent_reads.as_os_str(),
                OsStr::new("-o"),
                parent_extracted_fa.as_os_str(),
            ],
        )?;
    }

    if up_to_date(&[&parent_extracted_fa, &args.parents_ref], &[&parent_bam]) && !args.force {
        info!("Parent alignment already exists: {}", parent_bam.display());
    } else {
        run_pipe(
            &work,
            (
                "minimap2",
                vec![
                    OsStr::new("-t"),
                    OsStr::new(&args.threads_minimap2.to_string()),
                    OsStr::new("-ax"),
                    OsStr::new(args.preset.as_str()),
                    OsStr::new("--eqx"),
                    OsStr::new("--secondary=no"),
                    args.parents_ref.as_os_str(),
                    parent_extracted_fa.as_os_str(),
                ],
            ),
            (
                "samtools",
                vec![
                    OsStr::new("sort"),
                    OsStr::new("-@"),
                    OsStr::new(&args.threads_sort.to_string()),
                    OsStr::new("-o"),
                    parent_bam.as_os_str(),
                ],
            ),
        )?;
    }

    // 4) Regions present in F1 but NOT parent (control)
    if up_to_date(&[&f1_bam, &parent_bam], &[&regions_tsv]) && !args.force {
        info!("Regions TSV already exists: {}", regions_tsv.display());
    } else {
        run(
            &work,
            "klassify",
            &[
                OsStr::new("regions"),
                f1_bam.as_os_str(),
                parent_bam.as_os_str(),
            ],
        )?;
    }

    // 5) Refine the breakpoints (extract BAM segments overlapping regions, call breakpoints)
    let split_output = work.join("f1_classify.regions.split.fasta");
    if up_to_date(&[&regions_tsv, &f1_bam], &[&regions_fa, &split_output]) && !args.force {
        info!("Regions FASTA already exists: {}", regions_fa.display());
    } else {
        run(
            &work,
            "klassify",
            &[
                OsStr::new("extract-bam"),
                regions_tsv.as_os_str(),
                f1_bam.as_os_str(),
            ],
        )?;
        // klassify decides names internally; next call uses regions_fa
        run(
            &work,
            "klassify",
            &[
                OsStr::new("breakpoint"),
                kmers_bc.as_os_str(),
                regions_fa.as_os_str(),
            ],
        )?;
    }

    // Remap the split reads to parents to get a crisp ROI BAM
    if up_to_date(&[&split_output], &[&roi_bam]) && !args.force {
        info!("ROI BAM already exists: {}", roi_bam.display());
    } else {
        run_pipe(
            &work,
            (
                "minimap2",
                vec![
                    OsStr::new("-t"),
                    OsStr::new(&args.threads_minimap2.to_string()),
                    OsStr::new("-ax"),
                    OsStr::new(args.preset.as_str()),
                    OsStr::new("--eqx"),
                    OsStr::new("--secondary=no"),
                    args.parents_ref.as_os_str(),
                    split_output.as_os_str(),
                ],
            ),
            (
                "samtools",
                vec![
                    OsStr::new("sort"),
                    OsStr::new("-@"),
                    OsStr::new(&args.threads_sort.to_string()),
                    OsStr::new("-o"),
                    roi_bam.as_os_str(),
                ],
            ),
        )?;
    }

    // Cluster paired regions (capture stdout into ROI TSV)
    if up_to_date(&[&roi_bam], &[&roi_tsv]) && !args.force {
        info!("Paired regions TSV already exists: {}", roi_tsv.display());
    } else {
        run_to_file(
            &work,
            "klassify",
            &[OsStr::new("cluster-pairs"), roi_bam.as_os_str()],
            &roi_tsv,
        )?;
    }

    info!("Done.");
    info!(
        "Key outputs:\n  k-mers:          {}\n  F1 BAM:          {}\n  Parent BAM:      {}\n  Regions TSV:     {}\n  Regions FASTA:   {}\n  ROI BAM:         {}\n  ROI TSV:         {}\n  Paired Regions:  {}",
        kmers_bc.display(),
        f1_bam.display(),
        parent_bam.display(),
        regions_tsv.display(),
        regions_fa.display(),
        roi_bam.display(),
        roi_tsv.display(),
        paired_regions.display()
    );

    Ok(())
}

// ---------- helpers ----------

fn mkdir(p: &Path) -> Result<()> {
    if !p.exists() {
        fs::create_dir_all(p).with_context(|| format!("Failed to create {}", p.display()))?;
    }
    Ok(())
}

fn check_tool(tool: &str) -> Result<()> {
    match which(tool) {
        Ok(path) => {
            debug!("Found tool {} at {}", tool, path.display());
            Ok(())
        }
        Err(_) => bail!("Required tool not found in PATH: {}", tool),
    }
}

/// Return true if *all outputs* exist and are newer or equal to *all inputs*.
fn up_to_date(inputs: &[&Path], outputs: &[&Path]) -> bool {
    if outputs.is_empty() {
        return false;
    }
    let mut in_times = Vec::new();
    for i in inputs {
        match fs::metadata(i).and_then(|m| m.modified()) {
            Ok(t) => in_times.push(t),
            Err(_) => return false, // missing input? treat as not up-to-date to be safe
        }
    }
    let mut out_times = Vec::new();
    for o in outputs {
        match fs::metadata(o).and_then(|m| m.modified()) {
            Ok(t) => out_times.push(t),
            Err(_) => return false, // missing output -> not up to date
        }
    }
    if in_times.is_empty() || out_times.is_empty() {
        return false;
    }
    let max_in = in_times.into_iter().max().unwrap_or(SystemTime::UNIX_EPOCH);
    let min_out = out_times
        .into_iter()
        .min()
        .unwrap_or(SystemTime::UNIX_EPOCH);
    min_out >= max_in
}

fn run(cwd: &Path, program: &str, args: &[&OsStr]) -> Result<()> {
    info!("$ {} {}", program, join_args(args));
    let status = Command::new(program)
        .args(args)
        .current_dir(cwd)
        .status()
        .with_context(|| format!("Failed to spawn {}", program))?;
    if !status.success() {
        bail!("Command failed (exit {}): {}", status, program);
    }
    Ok(())
}

fn run_pipe(cwd: &Path, cmd1: (&str, Vec<&OsStr>), cmd2: (&str, Vec<&OsStr>)) -> Result<()> {
    info!(
        "$ {} {} | {} {}",
        cmd1.0,
        join_args(&cmd1.1),
        cmd2.0,
        join_args(&cmd2.1)
    );
    // first: spawn cmd1 with stdout piped
    let mut p1 = Command::new(cmd1.0)
        .args(&cmd1.1)
        .current_dir(cwd)
        .stdout(Stdio::piped())
        .spawn()
        .with_context(|| format!("Failed to spawn {}", cmd1.0))?;

    // second: spawn cmd2 with stdin from p1.stdout
    let stdout1 = p1
        .stdout
        .take()
        .ok_or_else(|| anyhow!("Failed to take stdout from {}", cmd1.0))?;
    let p2 = Command::new(cmd2.0)
        .args(&cmd2.1)
        .current_dir(cwd)
        .stdin(stdout1)
        .status()
        .with_context(|| format!("Failed to run {}", cmd2.0))?;

    // wait for p1 to complete
    let status1 = p1
        .wait()
        .with_context(|| format!("Failed to wait for {}", cmd1.0))?;

    if !status1.success() {
        bail!("Upstream command failed (exit {}): {}", status1, cmd1.0);
    }
    if !p2.success() {
        bail!("Downstream command failed (exit {}): {}", p2, cmd2.0);
    }
    Ok(())
}

fn run_to_file(cwd: &Path, program: &str, args: &[&OsStr], out_path: &Path) -> Result<()> {
    info!("$ {} {} > {}", program, join_args(args), out_path.display());
    let mut child = Command::new(program)
        .args(args)
        .current_dir(cwd)
        .stdout(Stdio::piped())
        .spawn()
        .with_context(|| format!("Failed to spawn {}", program))?;

    let mut reader = child
        .stdout
        .take()
        .ok_or_else(|| anyhow!("Failed to capture stdout from {}", program))?;

    let mut writer = BufWriter::new(
        File::create(out_path).with_context(|| format!("Cannot create {}", out_path.display()))?,
    );

    std::io::copy(&mut reader, &mut writer)
        .with_context(|| format!("Failed to write {}", out_path.display()))?;
    writer.flush().ok();

    let status = child
        .wait()
        .with_context(|| format!("Failed to wait for {}", program))?;
    if !status.success() {
        bail!("Command failed (exit {}): {}", status, program);
    }
    Ok(())
}

fn join_args(args: &[&OsStr]) -> String {
    args.iter()
        .map(|a| a.to_string_lossy().into_owned())
        .collect::<Vec<_>>()
        .join(" ")
}
