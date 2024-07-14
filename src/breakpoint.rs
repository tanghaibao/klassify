use crate::info::{load_kmer_db, map_kmer_to_file};
use crate::models::SingletonKmers;

use clap::Parser;
use log;
use needletail::{parse_fastx_file, Sequence};
use rayon::prelude::*;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

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
    // Classify the reads
    let mut reader = parse_fastx_file(fasta_file).expect("valid FASTA file");
    let file_prefix: &str = Path::new(fasta_file).file_name().unwrap().to_str().unwrap();
    let output_file = file_prefix.to_string() + ".classifications.bed";
    let mut writer = BufWriter::new(File::create(&output_file).unwrap());
    log::info!("Parsing reference");

    // Iterate through the reads
    let kmer_size = singleton_kmers.kmer_size;
    while let Some(record) = reader.next() {
        let record = record.expect("valid record");
        let seq = record.normalize(false);
        // Get the first part of the ID
        let id = String::from_utf8(record.id().to_vec())
            .unwrap()
            .split_whitespace()
            .next()
            .unwrap()
            .to_string();
        for (pos, kmer, _) in seq.bit_kmers(kmer_size, true) {
            if let Some(&file_index) = kmer_to_file.get(&kmer.0) {
                let to_write: String = format!(
                    "{}\t{}\t{}\t{}:{}",
                    id,
                    pos,
                    pos + kmer_size as usize,
                    singleton_kmers.fasta_files[file_index],
                    kmer.0
                );
                writeln!(writer, "{}", to_write).unwrap();
            }
        }
    }
    log::info!("Classifications written to `{}`", output_file);
}
