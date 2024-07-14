use crate::models::{ClassifyResults, SingletonKmers};
use bincode::deserialize_from;
use clap::Parser;
use log;
use needletail::{parse_fastx_file, Sequence};
use std::{
    collections::HashMap,
    fs::File,
    io::{BufReader, BufWriter, Write},
};

#[derive(Parser, Debug)]
pub struct ClassifyArgs {
    /// Bincode file
    pub bincode_file: String,
    /// Read file to classify
    pub reads_file: String,
}

pub fn classify(bincode_file: &str, reads_file: &str) {
    let reader = BufReader::new(File::open(bincode_file).unwrap());
    let singleton_kmers: SingletonKmers = deserialize_from(reader).unwrap();
    let kmer_size = singleton_kmers.kmer_size;
    log::info!(
        "Loaded singleton kmers (K={}) from `{}`",
        kmer_size,
        bincode_file
    );

    // Convert to kmer => file index
    let mut kmer_to_file = HashMap::new();
    for (file_index, kmer_set) in singleton_kmers.kmers.iter().enumerate() {
        for &kmer in kmer_set.iter() {
            kmer_to_file.insert(kmer, file_index);
        }
    }
    log::info!("Mapped kmers to files");
    // Classify the reads
    let mut reader = parse_fastx_file(reads_file).expect("valid reads file");
    let output_file = "read_classifications.txt";
    let mut writer = BufWriter::new(File::create(output_file).unwrap());
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
    while let Some(record) = reader.next() {
        let record = record.expect("valid record");
        let seq = record.normalize(false);
        let mut counts = vec![0; singleton_kmers.n()];
        for (_, kmer, _) in seq.bit_kmers(kmer_size, true) {
            if let Some(file_index) = kmer_to_file.get(&kmer.0) {
                counts[*file_index] += 1;
            }
        }
        let scaled_counts = counts
            .iter()
            .zip(singleton_kmers.scaling_factors.iter())
            .map(|(count, scaling_factor)| ((*count as f64) * scaling_factor) as i32)
            .collect::<Vec<_>>();
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
            scaled_counts,
        };
        let tag = results.tag(&singleton_kmers.fasta_files);
        let to_write = format!(
            "{}\t{}\t{}\t{}\t{}",
            results.id,
            results.seq_len,
            results.counts.iter().sum::<i32>(),
            tag,
            results
                .scaled_counts
                .iter()
                .map(|count| count.to_string())
                .collect::<Vec<_>>()
                .join("\t")
        );
        writeln!(writer, "{}", to_write).unwrap();
    }
    log::info!("Read classifications written to `{}`", output_file);
}
