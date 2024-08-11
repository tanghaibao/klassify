use log;
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::Path;

/// Discrete bin size to contract regions
pub const BINSIZE: u32 = 10_000;
/// Chain distance to merge regions
pub const CHAIN_DISTANCE: u32 = 2 * BINSIZE;
/// Flank size to extract from the region
pub const DEFAULT_FLANK_SIZE: i32 = BINSIZE as i32;
/// Maximum divergence
pub const MAX_DE: f32 = 0.01; // 1%

#[derive(Serialize, Deserialize)]
pub struct SingletonKmers {
    pub kmer_size: u8,
    pub fasta_files: Vec<String>,
    pub kmers: Vec<Vec<u64>>,
}

impl SingletonKmers {
    #[inline]
    pub fn n(&self) -> usize {
        self.fasta_files.len()
    }
}

pub struct ClassifyResults {
    pub id: String,
    pub seq_len: usize,
    pub counts: Vec<i32>,
}

impl ClassifyResults {
    pub fn tag(&self, fasta_files: &Vec<String>) -> String {
        let mut best_count = 0;
        let mut best_index = 0;
        let mut second_best_count = 0;
        let mut second_best_index = 0;
        let mut total = 0;
        for (i, &count) in self.counts.iter().enumerate() {
            total += count;
            if count > best_count {
                second_best_count = best_count;
                second_best_index = best_index;
                best_count = count;
                best_index = i;
            } else if count > second_best_count {
                second_best_count = count;
                second_best_index = i;
            }
        }
        let half = total / 2;
        if total == 0 {
            "Unclassified:0,0".to_string()
        } else if best_count + second_best_count > half {
            format!(
                "{},{}:{},{}",
                fasta_files[best_index],
                fasta_files[second_best_index],
                best_count * 100 / total,
                second_best_count * 100 / total
            )
        } else {
            format!(
                "Unclassified:{},{}",
                best_count * 100 / total,
                second_best_count * 100 / total
            )
        }
    }
}

/// Get basename
pub fn prefix(file_path: &str) -> String {
    Path::new(file_path)
        .file_name()
        .unwrap()
        .to_str()
        .unwrap()
        .to_string()
}

/// Get basename and up to the first dot in path
pub fn prefix_until_dot(file_path: &str) -> String {
    prefix(file_path).split('.').next().unwrap().to_string()
}

/// Run shell command
pub fn sh(command: &str) -> bool {
    log::info!("{}", command);
    let status = std::process::Command::new("sh")
        .arg("-c")
        .arg(command)
        .status()
        .expect("failed to execute process");
    status.success()
}

/// Determine if file a is newer than file b
fn is_newer_file(a: &str, b: &str) -> bool {
    let a_modified = fs::metadata(a).and_then(|m| m.modified()).ok();
    let b_modified = fs::metadata(b).and_then(|m| m.modified()).ok();

    match (a_modified, b_modified) {
        (Some(a_time), Some(b_time)) => a_time > b_time,
        _ => false,
    }
}

/// Check if any file in list a is newer than file in list b
pub fn need_update(a: Vec<String>, b: Vec<String>, warn: bool) -> bool {
    let should_update = b.iter().any(|x| !Path::new(x).exists())
        || b.iter()
            .all(|x| fs::metadata(x).map(|m| m.len() == 0).unwrap_or(false))
        || a.iter().any(|x| b.iter().any(|y| is_newer_file(x, y)));

    if !should_update && warn {
        log::info!("File `{}` found. Computation skipped.", b.join(", "));
    }

    should_update
}
