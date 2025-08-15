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
