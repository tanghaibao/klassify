use crate::info::{load_kmer_db, map_kmer_to_file};
use crate::models::{prefix, ClassifyResults, SingletonKmers};

use clap::Parser;
use log;
use needletail::{parse_fastx_file, Sequence};
use rayon::prelude::*;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
};

#[derive(Parser, Debug)]
pub struct ClassifyArgs {
    /// Bincode file
    pub bincode_file: String,
    /// Read file to classify
    pub reads_file: Vec<String>,
}

/// Classify reads based on unique (singleton) kmers.
pub fn classify(bincode_file: &str, reads_files: &Vec<String>) {
    let singleton_kmers = load_kmer_db(bincode_file);
    let kmer_to_file = map_kmer_to_file(&singleton_kmers);
    reads_files.par_iter().for_each(|reads_file| {
        classify_one(&singleton_kmers, &kmer_to_file, reads_file);
    });
}

fn classify_one(
    singleton_kmers: &SingletonKmers,
    kmer_to_file: &HashMap<u64, usize>,
    reads_file: &str,
) {
    // Classify the reads
    let mut reader = parse_fastx_file(reads_file).expect("valid reads file");
    let file_prefix = prefix(reads_file);
    let output_file = file_prefix + ".read_classifications.tsv";
    let mut writer = BufWriter::new(File::create(&output_file).unwrap());
    log::info!("Classifying reads");
    writeln!(
        writer,
        "ID\tLength\tKmers\tClassification\t{}",
        singleton_kmers
            .fasta_files
            .iter()
            .map(|count| count.to_string())
            .collect::<Vec<_>>()
            .join("\t")
    )
    .unwrap();

    // Iterate through the reads
    let kmer_size = singleton_kmers.kmer_size;
    while let Some(record) = reader.next() {
        let record = record.expect("valid record");
        let seq = record.normalize(false);
        let mut counts = vec![0; singleton_kmers.n()];
        for (_, kmer, _) in seq.bit_kmers(kmer_size, true) {
            if let Some(file_index) = kmer_to_file.get(&kmer.0) {
                counts[*file_index] += 1;
            }
        }
        // Get the first part of the ID
        let id = String::from_utf8(record.id().to_vec())
            .unwrap()
            .split_whitespace()
            .next()
            .unwrap()
            .to_string();
        let results = ClassifyResults {
            id,
            seq_len: record.seq().len(),
            counts,
        };
        let tag = results.tag(&singleton_kmers.fasta_files);
        let to_write = format!(
            "{}\t{}\t{}\t{}\t{}",
            results.id,
            results.seq_len,
            results.counts.iter().sum::<i32>(),
            tag,
            results
                .counts
                .iter()
                .map(|count| count.to_string())
                .collect::<Vec<_>>()
                .join("\t")
        );
        writeln!(writer, "{}", to_write).unwrap();
    }
    log::info!("Read classifications written to `{}`", output_file);
}
