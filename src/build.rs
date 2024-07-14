use crate::models::SingletonKmers;
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
    // Serialize the singleton kmers to a file
    let counts = singletons
        .iter()
        .map(|x| x.len() as f64)
        .collect::<Vec<_>>();
    let max_count = counts.iter().cloned().fold(0. / 0., f64::max);
    let scaling_factors = counts.iter().map(|&x| max_count / x).collect();
    let singleton_kmers = SingletonKmers {
        kmer_size: KMER_SIZE,
        fasta_files: fasta_files.clone(),
        kmers: singletons,
        scaling_factors,
    };
    let output_file = "singleton_kmers.bc";
    let writer = BufWriter::new(File::create(output_file).unwrap());
    serialize_into(writer, &singleton_kmers).expect("serialization to succeed");
    log::info!("Singleton kmers written to `{}`", output_file);
}

fn get_kmers(fasta_file: &str) -> HashSet<u64> {
    let mut reader = parse_fastx_file(fasta_file).expect("valid FASTA file");
    let mut kmer_set = HashSet::new();
    while let Some(record) = reader.next() {
        let record = record.expect("valid record");
        let seq = record.normalize(false);
        for (_, kmer, _) in seq.bit_kmers(KMER_SIZE, true) {
            kmer_set.insert(kmer.0);
        }
    }
    log::info!("{}: {} kmers found", fasta_file, kmer_set.len());
    kmer_set
}
