use bincode::{Decode, Encode};
use log::info;
use rust_htslib::bam;
use std::fs::metadata;
use std::path::Path;

/// Discrete bin size to contract regions
pub const BINSIZE: u32 = 10_000;

/// Chain distance to merge regions
pub const CHAIN_DISTANCE: u32 = 2 * BINSIZE;

/// Flank size to extract from the region
pub const DEFAULT_FLANK_SIZE: i32 = BINSIZE as i32;

/// Maximum divergence
pub const MAX_DE: f32 = 0.01; // 1%

/// Number of threads used to read BAM files
const N_THREADS_READ_BAM: usize = 32;

#[derive(Decode, Encode)]
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
    pub fn tag(&self, fasta_files: &[String]) -> String {
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

/// Get basename and up to the first dot in the path
pub fn prefix_until_dot(file_path: &str) -> String {
    prefix(file_path).split('.').next().unwrap().to_string()
}

/// Run shell command
pub fn sh(command: &str) -> bool {
    info!("{}", command);
    let status = std::process::Command::new("sh")
        .arg("-c")
        .arg(command)
        .status()
        .expect("failed to execute process");
    status.success()
}

/// Determine if file `a` is newer than the file `b`
fn is_newer_file<S: AsRef<Path>, T: AsRef<Path>>(a: S, b: T) -> bool {
    let a_modified = metadata(a).and_then(|m| m.modified()).ok();
    let b_modified = metadata(b).and_then(|m| m.modified()).ok();

    match (a_modified, b_modified) {
        (Some(a_time), Some(b_time)) => a_time > b_time,
        _ => false,
    }
}

/// Check if any file in the list `a` is newer than any file in the list `b`
pub fn need_update<S: AsRef<Path>, T: AsRef<Path>>(a: &[S], b: &[T], warn: bool) -> bool {
    let should_update = b.iter().any(|x| !Path::new(x.as_ref()).exists())
        || b.iter()
            .all(|x| metadata(x).map(|m| m.len() == 0).unwrap_or(false))
        || a.iter().any(|x| b.iter().any(|y| is_newer_file(x, y)));

    if !should_update && warn {
        let files: Vec<String> = b
            .iter()
            .map(|x| format!("{}", x.as_ref().display()))
            .collect();
        info!("File `{}` found. Computation skipped.", files.join(", "));
    }

    should_update
}

/// Build BAM index if needed
pub fn index_bam(bam_file: &str) {
    let bam_index = bam_file.to_string() + ".bai";
    if !need_update(&[bam_file], &[&bam_index], false) {
        info!("BAM index `{bam_index}` already exists. Skipping indexing.");
        return;
    }
    let n_threads = N_THREADS_READ_BAM;
    info!("Indexing `{bam_file}` (n_threads={n_threads})");
    _ = bam::index::build(bam_file, None, bam::index::Type::Bai, n_threads as u32)
}

#[cfg(test)]
mod tests {
    use super::*;
    use filetime::FileTime;
    use std::fs::{write, File};
    use std::time::{Duration, SystemTime};
    use tempfile::TempDir;

    #[test]
    fn prefix_extracts_basename() {
        assert_eq!(prefix("/tmp/foo/bar.baz"), "bar.baz");
    }

    #[test]
    fn prefix_until_dot_truncates_correctly() {
        assert_eq!(prefix_until_dot("/tmp/foo/bar.baz"), "bar");
        assert_eq!(prefix_until_dot("no_dot"), "no_dot"); // no dot
    }

    #[test]
    fn sh_success_returns_true() {
        assert!(sh("true"));
    }

    #[test]
    fn sh_failure_returns_false() {
        assert!(!sh("exit 1"));
    }

    fn touch(path: &Path, offset_secs: i64) {
        // Helper to set mtime relative to now
        let now = SystemTime::now()
            .checked_add(Duration::from_secs(offset_secs.unsigned_abs()))
            .unwrap();
        let ft = FileTime::from_system_time(now);
        filetime::set_file_times(path, ft, ft).unwrap();
    }

    #[test]
    fn need_update_when_target_missing() {
        let dir = TempDir::new().unwrap();
        let src = dir.path().join("src");
        File::create(&src).unwrap();
        let missing = dir.path().join("missing");
        assert!(need_update(&[&src], &[&missing], false));
    }

    #[test]
    fn need_update_when_target_is_zero_bytes() {
        let dir = TempDir::new().unwrap();
        let src = dir.path().join("src");
        let tgt = dir.path().join("tgt");
        write(&src, b"data").unwrap();
        File::create(&tgt).unwrap(); // zero length
        assert!(need_update(&[&src], &[tgt], false));
    }

    #[test]
    fn need_update_when_source_newer_than_target() {
        let dir = TempDir::new().unwrap();
        let src = dir.path().join("src");
        let tgt = dir.path().join("tgt");
        write(&src, b"old").unwrap();
        touch(&src, -10); // 10 s ago
        write(&tgt, b"new").unwrap(); // newer
        assert!(need_update(&[&src], &[tgt], false));
    }

    #[test]
    fn need_update_false_when_up_to_date() {
        let dir = TempDir::new().unwrap();
        let src = dir.path().join("src");
        let tgt = dir.path().join("tgt");
        write(&tgt, b"old").unwrap();
        touch(&tgt, -10); // 10 s ago
        write(&src, b"new").unwrap(); // newer
        assert!(!need_update(&[&src], &[tgt], false));
    }
}
