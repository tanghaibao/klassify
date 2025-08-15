use crate::utils::SingletonKmers;
use bincode::{config, encode_into_std_write};
use clap::Parser;
use log::info;
use needletail::{parse_fastx_file, Sequence};
use rayon::prelude::*;
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::BufWriter,
};

/// K-mer length
const KMER_SIZE: u8 = 24;
/// File name for singleton kmers
const SINGLETON_KMERS: &str = "kmers.bc";

#[derive(Parser, Debug)]
pub struct BuildArgs {
    /// Input FASTA files
    pub fasta_files: Vec<String>,
    /// Output file
    #[clap(short, long, default_value = SINGLETON_KMERS)]
    pub output_file: String,
    /// K-mer size
    #[clap(short, long, default_value_t = KMER_SIZE)]
    pub kmer_size: u8,
}

/// Sequence records with an ID and sequence
struct SeqRecord {
    pub id: String,
    pub seq: Vec<u8>,
}

/// Convert FASTA files to singleton k-mers
pub fn build(fasta_files: &Vec<String>, output_file: &str, kmer_size: u8) {
    let all_records = fasta_files
        .par_iter()
        .flat_map(|fasta_file| get_sequences(fasta_file))
        .collect::<Vec<_>>();
    let all_sets = all_records
        .par_iter()
        .map(|seq_record| get_kmers(seq_record, kmer_size))
        .collect::<Vec<_>>();
    // Identify all the kmers that appear once and only once in all the files
    let mut kmer_counts = HashMap::new();
    for kmer_set in all_sets.iter() {
        for kmer in kmer_set.iter() {
            *kmer_counts.entry(*kmer).or_insert(0) += 1;
        }
    }
    info!("Total unique kmers: {}", kmer_counts.len());
    let singleton_kmers = kmer_counts
        .into_iter()
        .filter(|(_, count)| *count == 1)
        .map(|(kmer, _)| kmer)
        .collect::<HashSet<_>>();
    info!("Singleton kmers: {}", singleton_kmers.len());
    let ids = all_records
        .iter()
        .map(|seq_record| seq_record.id.clone())
        .collect::<Vec<_>>();
    // Find the unique kmers in each file
    let kmers = (all_records, all_sets)
        .into_par_iter()
        .map(|(seq_record, kmer_set)| {
            let singleton_kmers_per_file = kmer_set
                .intersection(&singleton_kmers)
                .cloned()
                .collect::<Vec<_>>();
            info!(
                "{}: {} singleton kmers found",
                seq_record.id,
                singleton_kmers_per_file.len()
            );
            singleton_kmers_per_file
        })
        .collect::<Vec<_>>();
    // Serialize the singleton kmers to a file
    let singleton_kmers = SingletonKmers {
        kmer_size,
        ids,
        kmers,
    };
    let mut writer = BufWriter::new(File::create(output_file).unwrap());
    encode_into_std_write(&singleton_kmers, &mut writer, config::standard())
        .expect("Failed to serialize");
    info!("Singleton kmers written to `{}`", output_file);
}

/// Get kmers from a FASTA file
fn get_sequences(fasta_file: &str) -> Vec<SeqRecord> {
    let mut reader = parse_fastx_file(fasta_file).expect("valid FASTA file");
    let mut records = Vec::new();
    while let Some(record) = reader.next() {
        let record = record.expect("valid record");
        let id = String::from_utf8(record.id().to_vec()).expect("valid UTF-8");
        let seq = record.normalize(false).to_vec();
        records.push(SeqRecord { id, seq });
    }
    info!("{}: {} records found", fasta_file, records.len());
    records
}

/// Get kmers from a FASTA file
fn get_kmers(seq_record: &SeqRecord, kmer_size: u8) -> HashSet<u64> {
    let mut kmer_set = HashSet::new();
    for (_, kmer, _) in seq_record.seq.bit_kmers(kmer_size, true) {
        kmer_set.insert(kmer.0);
    }
    info!("{}: {} kmers found", seq_record.id, kmer_set.len());
    kmer_set
}

#[cfg(test)]
mod tests {
    use super::*;
    use bincode::{config, decode_from_std_read};
    use std::fs::{self, File};
    use std::io::BufReader;
    use std::path::{Path, PathBuf};
    use tempfile::tempdir;

    /// Write a small FASTA file into `dir` and return its path.
    fn write_fasta(dir: &Path, name: &str, contents: &str) -> PathBuf {
        let p = dir.join(name);
        fs::write(&p, contents).expect("write fasta");
        p
    }

    #[test]
    fn test_get_kmers_canonicalization() {
        // Sequence contains both ACG and its RC (CGT); with k=3 and canonical=true,
        // these should collapse to a single k-mer in the set, and similarly CGC/GCG collapse.
        let rec = SeqRecord {
            id: "x".to_string(),
            seq: b"ACGCGT".to_vec(),
        };
        let set = get_kmers(&rec, 3);
        // 6bp sequence, 4 windows, but canonical pairs collapse => {ACG, CGC} = 2
        assert_eq!(
            set.len(),
            2,
            "canonicalized 3-mers should collapse RC pairs"
        );
    }

    #[test]
    fn test_build_finds_singletons_and_writes_bincode() {
        let tmp = tempdir().expect("tmpdir");
        let out = tmp.path().join("kmers.bc");

        // File 1 has two records; File 2 has one record.
        // With k=3 (canonical), the only global singleton should be TAA (present only in r2).
        let f1 = write_fasta(
            tmp.path(),
            "a.fa",
            r#">r1
ACGTAC
>r2
ACGTAA
"#,
        );
        let f2 = write_fasta(
            tmp.path(),
            "b.fa",
            r#">r3
ACGTAC
"#,
        );

        let fasta_files = vec![
            f1.to_string_lossy().to_string(),
            f2.to_string_lossy().to_string(),
        ];

        // Run the pipeline with k=3 to keep numbers tiny & deterministic.
        build(&fasta_files, out.to_str().unwrap(), 3);

        // The file should exist and be non-empty.
        assert!(out.exists(), "output bincode file should be written");
        let meta = fs::metadata(&out).expect("stat out");
        assert!(meta.len() > 0, "output file should be non-empty");

        // Decode and inspect structure/content.
        let mut reader = BufReader::new(File::open(&out).expect("open out"));
        let decoded: SingletonKmers =
            decode_from_std_read(&mut reader, config::standard()).expect("bincode decode");

        // Sanity checks.
        assert_eq!(decoded.kmer_size, 3, "kmer_size should be what we passed");
        assert_eq!(
            decoded.ids.len(),
            decoded.kmers.len(),
            "ids and kmers should be aligned"
        );
        assert_eq!(
            decoded.ids.len(),
            3,
            "we had 3 sequence records across the two files"
        );

        // We don't rely on record order because rayon's parallel iterator doesn't guarantee it.
        // Instead, confirm exactly one record has a singleton and it belongs to r2.
        let total_singletons: usize = decoded.kmers.iter().map(|v| v.len()).sum();
        assert_eq!(
            total_singletons, 1,
            "there should be exactly one singleton k-mer overall"
        );

        // Find r2 and assert it holds the singleton.
        let mut found_r2 = false;
        for (id, kmers) in decoded.ids.iter().zip(decoded.kmers.iter()) {
            if id == "r2" {
                found_r2 = true;
                assert_eq!(kmers.len(), 1, "r2 should have the only singleton (TAA)");
            } else {
                assert!(
                    kmers.is_empty(),
                    "records other than r2 should have zero singletons"
                );
            }
        }
        assert!(found_r2, "decoded ids should include r2");
    }

    #[test]
    fn test_buildargs_parse_defaults() {
        // Check that clap defaults match the constants.
        let args = BuildArgs::parse_from(["prog", "a.fa", "b.fa"]);
        assert_eq!(args.output_file, SINGLETON_KMERS);
        assert_eq!(args.kmer_size, KMER_SIZE);
        assert_eq!(args.fasta_files.len(), 2);
        assert!(args.fasta_files.contains(&"a.fa".to_string()));
        assert!(args.fasta_files.contains(&"b.fa".to_string()));
    }
}
