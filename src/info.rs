use crate::models::SingletonKmers;

use bincode::deserialize_from;
use clap::Parser;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;

#[derive(Parser, Debug)]
pub struct InfoArgs {
    /// Bincode file
    pub bincode_file: String,
}

/// Load the bincode file and print some information about it.
pub fn load_kmer_db(bincode_file: &str) -> SingletonKmers {
    let reader = BufReader::new(File::open(bincode_file).unwrap());
    let singleton_kmers: SingletonKmers = deserialize_from(reader).unwrap();
    log::info!(
        "Loaded singleton kmers (K={}) from `{}`",
        singleton_kmers.kmer_size,
        bincode_file
    );
    singleton_kmers
}

pub fn map_kmer_to_file(singleton_kmers: &SingletonKmers) -> HashMap<u64, usize> {
    // Convert to kmer => file index
    let mut kmer_to_file = HashMap::new();
    for (file_index, kmer_set) in singleton_kmers.kmers.iter().enumerate() {
        for &kmer in kmer_set.iter() {
            kmer_to_file.insert(kmer, file_index);
        }
    }
    log::info!("Mapped kmers to files");
    kmer_to_file
}

pub fn info(bincode_file: &str) {
    let singleton_kmers = load_kmer_db(bincode_file);
    println!("Kmer size: {}", singleton_kmers.kmer_size);
    println!("Number of fasta files: {}", singleton_kmers.n());
    println!(
        "Number of kmers: {}",
        singleton_kmers.kmers.iter().map(|x| x.len()).sum::<usize>()
    );
    for (i, fasta_file) in singleton_kmers.fasta_files.iter().enumerate() {
        println!(
            "  {}: {} ({} mers)",
            i + 1,
            fasta_file,
            singleton_kmers.kmers[i].len()
        );
    }
}
