use crate::models::{prefix, SingletonKmers};
use bincode::serialize_into;
use clap::Parser;
use needletail::{parse_fastx_file, Sequence};
use rayon::prelude::*;
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::BufWriter,
};

const KMER_SIZE: u8 = 24;
const SINGLETON_KMERS: &str = "singleton_kmers.bc";

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

/// Convert FASTA files to singleton k-mers
pub fn build(fasta_files: &Vec<String>, output_file: &str, kmer_size: u8) {
    let all_sets = fasta_files
        .par_iter()
        .map(|fasta_file| get_kmers(fasta_file, kmer_size))
        .collect::<Vec<_>>();
    // Identify all the kmers that appear once and only once in all the files
    let mut kmer_counts = HashMap::new();
    for kmer_set in all_sets.iter() {
        for kmer in kmer_set.iter() {
            *kmer_counts.entry(*kmer).or_insert(0) += 1;
        }
    }
    log::info!("Total unique kmers: {}", kmer_counts.len());
    let singleton_kmers = kmer_counts
        .into_iter()
        .filter(|(_, count)| *count == 1)
        .map(|(kmer, _)| kmer)
        .collect::<HashSet<_>>();
    log::info!("Singleton kmers: {}", singleton_kmers.len());
    // Find the unique kmers in each file
    let singletons = (fasta_files, all_sets)
        .into_par_iter()
        .map(|(fasta_file, kmer_set)| {
            let singleton_kmers_per_file = kmer_set
                .intersection(&singleton_kmers)
                .cloned()
                .collect::<Vec<_>>();
            log::info!(
                "{}: {} singleton kmers found",
                fasta_file,
                singleton_kmers_per_file.len()
            );
            singleton_kmers_per_file
        })
        .collect::<Vec<_>>();
    let fasta_files = fasta_files.iter().map(|x| prefix(x)).collect::<Vec<_>>();
    // Serialize the singleton kmers to a file
    let singleton_kmers = SingletonKmers {
        kmer_size,
        fasta_files,
        kmers: singletons,
    };
    let writer = BufWriter::new(File::create(output_file).unwrap());
    serialize_into(writer, &singleton_kmers).expect("serialization to succeed");
    log::info!("Singleton kmers written to `{}`", output_file);
}

/// Get kmers from a FASTA file
fn get_kmers(fasta_file: &str, kmer_size: u8) -> HashSet<u64> {
    let mut reader = parse_fastx_file(fasta_file).expect("valid FASTA file");
    let mut kmer_set = HashSet::new();
    while let Some(record) = reader.next() {
        let record = record.expect("valid record");
        let seq = record.normalize(false);
        for (_, kmer, _) in seq.bit_kmers(kmer_size, true) {
            kmer_set.insert(kmer.0);
        }
    }
    log::info!("{}: {} kmers found", fasta_file, kmer_set.len());
    kmer_set
}
