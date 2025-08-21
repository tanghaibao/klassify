use crate::tools::info::{load_kmer_db, map_kmer_to_file};
use crate::utils::{prefix_until_dot, SingletonKmers};

use clap::Parser;
use log::{info, warn};
use needletail::{parse_fastx_file, Sequence};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

// =============================
// CLI
// =============================

const DEFAULT_KMER_THRESHOLD: usize = 30;
const DEFAULT_REGION_SEPARATOR: &str = "@";

#[derive(Parser, Debug)]
pub struct BreakpointArgs {
    /// Bincode file (singleton k-mers)
    pub bincode_file: String,
    /// FASTA/FASTQ files to process (gz ok)
    pub fasta_files: Vec<String>,

    /// If set, do NOT emit BED; write split FASTA (*.split.fasta) instead
    #[clap(long)]
    pub split: bool,

    /// Minimum number of k-mers supporting each side of the breakpoint (split mode only)
    #[clap(long, default_value_t = DEFAULT_KMER_THRESHOLD)]
    pub kmer_threshold: usize,

    /// Region separator used in read headers, e.g. "A@B@read123" (split mode only)
    #[clap(long, default_value = DEFAULT_REGION_SEPARATOR)]
    pub region_sep: String,
}

/// Public entry (call this from your main)
pub fn breakpoint(args: &BreakpointArgs) {
    let singleton_kmers = load_kmer_db(&args.bincode_file);
    let kmer_to_file = map_kmer_to_file(&singleton_kmers);

    if args.split {
        // Precompute accession prefixes once
        let accn_prefix: Vec<String> = singleton_kmers
            .ids
            .iter()
            .map(|s| prefix_until_dot(s))
            .collect();

        args.fasta_files.par_iter().for_each(|f| {
            split_reads_one(
                &singleton_kmers,
                &kmer_to_file,
                &accn_prefix,
                f,
                args.kmer_threshold,
                args.region_sep.chars().next().unwrap_or('@'),
            );
        });
    } else {
        // Original BED generation path (unchanged)
        args.fasta_files.par_iter().for_each(|f| {
            breakpoint_one_bed(&singleton_kmers, &kmer_to_file, f);
        });
    }
}

// =============================
// Original BED path (as before)
// =============================

fn breakpoint_one_bed(
    singleton_kmers: &SingletonKmers,
    kmer_to_file: &HashMap<u64, usize>,
    fasta_file: &str,
) {
    let mut reader = parse_fastx_file(fasta_file).expect("valid FASTA/FASTQ");
    let output_file = Path::new(fasta_file).with_extension("bed");
    let mut writer = BufWriter::new(File::create(&output_file).unwrap());
    info!("Parsing `{}` → `{}`", fasta_file, output_file.display());

    let kmer_size = singleton_kmers.kmer_size;
    while let Some(record) = reader.next() {
        let record = record.expect("valid record");
        let seq = record.normalize(false);

        // Get the first token of the ID
        let id = String::from_utf8(record.id().to_vec())
            .unwrap_or_else(|_| "unknown".into())
            .split_whitespace()
            .next()
            .unwrap_or("unknown")
            .to_string();

        for (pos, kmer, _) in seq.bit_kmers(kmer_size, true) {
            if let Some(&file_index) = kmer_to_file.get(&kmer.0) {
                let to_write: String = format!(
                    "{}\t{}\t{}\t{}:{}",
                    id,
                    pos,
                    pos + kmer_size as usize,
                    prefix_until_dot(&singleton_kmers.ids[file_index]),
                    kmer.0
                );
                writeln!(writer, "{to_write}").unwrap();
            }
        }
    }
    info!("BED written to `{}`", output_file.display());
}

// =============================
// Split mode (Python logic port)
// =============================

#[derive(Clone, Copy)]
struct Hit {
    start: usize,
    end: usize,
    file_idx: usize, // index into singleton_kmers.ids
}

fn split_reads_one(
    singleton_kmers: &SingletonKmers,
    kmer_to_file: &HashMap<u64, usize>,
    accn_prefix: &[String],
    fasta_file: &str,
    kmer_threshold: usize,
    region_sep: char,
) {
    let mut reader = parse_fastx_file(fasta_file).expect("valid FASTA/FASTQ");
    let output_file = Path::new(fasta_file).with_extension("split.fasta");
    let mut writer = BufWriter::new(File::create(&output_file).expect("create output"));
    let k = singleton_kmers.kmer_size as usize;

    info!(
        "Processing (split) `{}` -> `{}`",
        fasta_file,
        output_file.display()
    );

    let mut n_success = 0usize;
    let mut n_fail = 0usize;
    let mut n_not_enough = 0usize;

    while let Some(rec) = reader.next() {
        let rec = rec.expect("valid record");
        let id = String::from_utf8(rec.id().to_vec()).unwrap_or_default();

        // Expect "ra@rb@orig_read"
        let mut parts = id.splitn(3, region_sep);
        let ra = parts.next().unwrap_or("");
        let rb = parts.next().unwrap_or("");
        let orig_read = parts.next().unwrap_or(&id);

        if ra.is_empty() || rb.is_empty() {
            warn!(
                "Read `{}` missing `ra{sep}rb{sep}orig` format; skipped",
                id,
                sep = region_sep
            );
            n_fail += 1;
            continue;
        }

        // Collect k-mer hits
        let norm = rec.normalize(false);
        let mut hits: Vec<Hit> = Vec::new();
        for (pos, bk, _) in norm.bit_kmers(singleton_kmers.kmer_size, true) {
            if let Some(&fi) = kmer_to_file.get(&bk.0) {
                hits.push(Hit {
                    start: pos,
                    end: pos + k,
                    file_idx: fi,
                });
            }
        }
        if hits.len() < 2 {
            n_fail += 1;
            continue;
        }

        match choose_breakpoint(&hits, accn_prefix, ra, rb, kmer_threshold) {
            BreakResult::Ok {
                idx,
                left_label,
                right_label,
                ..
            } => {
                if idx + 1 >= hits.len() {
                    n_fail += 1;
                    continue;
                }
                let left_end = hits[idx].end;
                let right_start = hits[idx + 1].start;
                if left_end >= right_start {
                    n_fail += 1;
                    continue;
                }

                let seq = rec.seq();
                let end = seq.len();
                let mid = (left_end + right_start) / 2;
                if mid == 0 || mid >= end {
                    n_fail += 1;
                    continue;
                }

                let left_id = format!("{}|{}|0-{}", orig_read, left_label, mid);
                let right_id = format!("{}|{}|{}-{}", orig_read, right_label, mid, end);
                write_fasta(&mut writer, &left_id, &seq[..mid]);
                write_fasta(&mut writer, &right_id, &seq[mid..]);
                n_success += 1;
            }
            BreakResult::NotEnoughKmers => {
                n_not_enough += 1;
            }
            BreakResult::Fail => {
                n_fail += 1;
            }
        }
    }

    info!(
        "Summary (split) `{}`: SUCCESS={} NOT_ENOUGH_KMERS={} FAIL={}",
        fasta_file, n_success, n_not_enough, n_fail
    );
}

enum BreakResult<'a> {
    Ok {
        idx: usize,
        left_label: &'a str,
        right_label: &'a str,
    },
    NotEnoughKmers,
    Fail,
}

fn choose_breakpoint<'a>(
    hits: &[Hit],
    accn_prefix: &[String],
    ra: &'a str,
    rb: &'a str,
    kmer_threshold: usize,
) -> BreakResult<'a> {
    let n = hits.len();
    if n == 0 {
        return BreakResult::Fail;
    }
    let mut eq_ra = vec![false; n];
    let mut eq_rb = vec![false; n];
    for (i, h) in hits.iter().enumerate() {
        let acc = &accn_prefix[h.file_idx];
        eq_ra[i] = acc == ra;
        eq_rb[i] = acc == rb;
    }

    // prefix counts
    let mut pref_a = vec![0usize; n];
    let mut pref_b = vec![0usize; n];
    pref_a[0] = eq_ra[0] as usize;
    pref_b[0] = eq_rb[0] as usize;
    for i in 1..n {
        pref_a[i] = pref_a[i - 1] + eq_ra[i] as usize;
        pref_b[i] = pref_b[i - 1] + eq_rb[i] as usize;
    }

    // suffix counts (after i)
    let mut suf_a = vec![0usize; n];
    let mut suf_b = vec![0usize; n];
    for i in (1..n).rev() {
        suf_a[i - 1] = suf_a[i] + eq_ra[i] as usize;
        suf_b[i - 1] = suf_b[i] + eq_rb[i] as usize;
    }

    // ab = pref_a + suf_b; ba = pref_b + suf_a
    let mut ab_max = 0usize;
    let mut ab_idx = 0usize;
    let mut ba_max = 0usize;
    let mut ba_idx = 0usize;
    for i in 0..n {
        let ab = pref_a[i] + suf_b[i];
        if ab > ab_max {
            ab_max = ab;
            ab_idx = i;
        }
        let ba = pref_b[i] + suf_a[i];
        if ba > ba_max {
            ba_max = ba;
            ba_idx = i;
        }
    }

    // choose orientation
    let (idx, left_label, right_label, count_left, count_right) = if ab_max >= ba_max {
        (ab_idx, ra, rb, pref_a[ab_idx], suf_b[ab_idx])
    } else {
        (ba_idx, rb, ra, pref_b[ba_idx], suf_a[ba_idx])
    };

    if idx + 1 >= n {
        return BreakResult::Fail; // at end
    }
    if count_left < kmer_threshold || count_right < kmer_threshold {
        return BreakResult::NotEnoughKmers;
    }
    BreakResult::Ok {
        idx,
        left_label,
        right_label,
    }
}

fn write_fasta<W: Write>(writer: &mut W, id: &str, seq: &[u8]) {
    // uppercase without per-base allocation in fmt
    let mut up = Vec::with_capacity(seq.len());
    up.extend(seq.iter().map(|b| b.to_ascii_uppercase()));
    writeln!(writer, ">{}", id).unwrap();
    writer.write_all(&up).unwrap();
    writeln!(writer).unwrap();
}

// =============================
// Tests
// =============================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn breakpoint_pref_suf_orientation_ab() {
        // Map file_idx -> labels: 0 => "A", 1 => "B"
        let accn_prefix = vec!["A".to_string(), "B".to_string()];

        // hits: AAAA BBBB (break should be at i=3)
        let hits: Vec<Hit> = (0..4)
            .map(|i| Hit {
                start: i,
                end: i + 1,
                file_idx: 0,
            })
            .chain((4..8).map(|i| Hit {
                start: i,
                end: i + 1,
                file_idx: 1,
            }))
            .collect();

        match choose_breakpoint(&hits, &accn_prefix, "A", "B", 1) {
            BreakResult::Ok {
                idx,
                left_label,
                right_label,
                ..
            } => {
                assert_eq!(idx, 3);
                assert_eq!(left_label, "A");
                assert_eq!(right_label, "B");
            }
            _ => panic!("expected Ok"),
        }
    }

    #[test]
    fn breakpoint_orientation_ba() {
        // BBBB AAAA -> orientation should flip (left=B, right=A)
        let accn_prefix = vec!["A".to_string(), "B".to_string()];
        let hits: Vec<Hit> = (0..4)
            .map(|i| Hit {
                start: i,
                end: i + 1,
                file_idx: 1,
            })
            .chain((4..8).map(|i| Hit {
                start: i,
                end: i + 1,
                file_idx: 0,
            }))
            .collect();

        match choose_breakpoint(&hits, &accn_prefix, "A", "B", 1) {
            BreakResult::Ok {
                idx,
                left_label,
                right_label,
                ..
            } => {
                assert_eq!(idx, 3);
                assert_eq!(left_label, "B");
                assert_eq!(right_label, "A");
            }
            _ => panic!("expected Ok"),
        }
    }

    #[test]
    fn breakpoint_not_enough_kmers() {
        let accn_prefix = vec!["A".to_string(), "B".to_string()];
        let hits: Vec<Hit> = (0..2)
            .map(|i| Hit {
                start: i,
                end: i + 1,
                file_idx: 0,
            })
            .chain((2..4).map(|i| Hit {
                start: i,
                end: i + 1,
                file_idx: 1,
            }))
            .collect();
        match choose_breakpoint(&hits, &accn_prefix, "A", "B", 10) {
            BreakResult::NotEnoughKmers => {}
            _ => panic!("expected NotEnoughKmers"),
        }
    }

    #[test]
    fn breakpoint_fail_at_end() {
        // All A => best idx will be at end (no idx+1), should Fail
        let accn_prefix = vec!["A".to_string()];
        let hits: Vec<Hit> = (0..5)
            .map(|i| Hit {
                start: i,
                end: i + 1,
                file_idx: 0,
            })
            .collect();
        match choose_breakpoint(&hits, &accn_prefix, "A", "B", 1) {
            BreakResult::Fail => {}
            _ => panic!("expected Fail"),
        }
    }
}
