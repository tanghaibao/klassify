use clap::Parser;
use needletail::{parse_fastx_file, Sequence};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

const KMER_SIZE: u8 = 24;

#[derive(Parser, Debug)]
pub struct BuildArgs {
    /// Input FASTA files
    pub fasta_files: Vec<String>,
}

pub fn build(fasta_files: &Vec<String>) {
    let all_sets = fasta_files
        .par_iter()
        .map(|fasta_file| get_kmers(fasta_file))
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
    let _ = (fasta_files, all_sets)
        .into_par_iter()
        .map(|(fasta_file, kmer_set)| {
            let singleton_kmers_per_file = kmer_set
                .intersection(&singleton_kmers)
                .cloned()
                .collect::<HashSet<_>>();
            log::info!(
                "{}: {} singleton kmers found",
                fasta_file,
                singleton_kmers_per_file.len()
            );
            singleton_kmers_per_file
        })
        .collect::<Vec<_>>();
}

fn get_kmers(fasta_file: &str) -> HashSet<u64> {
    let mut reader = parse_fastx_file(fasta_file).expect("valid FASTA file");
    let mut kmer_set = HashSet::new();
    while let Some(record) = reader.next() {
        let record = record.expect("valid record");
        let seq = record.normalize(false);
        for (_, kmer, _) in seq.bit_kmers(KMER_SIZE, true) {
            kmer_set.insert(kmer.0);
            // let kmer_seq = bitmer_to_bytes(kmer);
            // println!(
            //     "{} {:?}",
            //     String::from_utf8(kmer_seq).expect("valid UTF-8"),
            //     kmer.0,
            // );
        }
    }
    log::info!("{}: {} kmers found", fasta_file, kmer_set.len());
    kmer_set
}
