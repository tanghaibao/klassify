use crate::info::{load_kmer_db, map_kmer_to_file};
use crate::models::SingletonKmers;

use clap::Parser;
use rayon::prelude::*;
use std::collections::HashMap;

#[derive(Parser, Debug)]
pub struct BreakpointArgs {
    /// Bincode file
    pub bincode_file: String,
    /// FASTA files to detect breakpoint
    pub fasta_files: Vec<String>,
}

pub fn breakpoint(bincode_file: &str, fasta_files: &Vec<String>) {
    let singleton_kmers = load_kmer_db(bincode_file);
    let kmer_to_file = map_kmer_to_file(&singleton_kmers);
    fasta_files.par_iter().for_each(|fasta_file| {
        breakpoint_one(&singleton_kmers, &kmer_to_file, fasta_file);
    });
}

fn breakpoint_one(
    singleton_kmers: &SingletonKmers,
    kmer_to_file: &HashMap<u64, usize>,
    fasta_file: &str,
) {
}
