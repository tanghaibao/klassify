use crate::info::{load_kmer_db, map_kmer_to_file};
use crate::models::{prefix, ClassifyResults, SingletonKmers};

use clap::Parser;
use log;
use needletail::{parse_fastx_file, Sequence};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::{collections::HashMap, io::BufWriter};

#[derive(Parser, Debug)]
pub struct ClassifyArgs {
    /// Bincode file
    pub bincode_file: String,
    /// Read file to classify
    pub reads_file: Vec<String>,
    /// Output directory
    #[clap(short, long)]
    pub output_dir: String,
    /// Prefix length
    #[clap(short, long, default_value_t = 7)]
    pub prefix_length: usize,
}

type ReadClassification = Vec<String>;

/// Classify reads based on unique (singleton) kmers.
pub fn classify(
    bincode_file: &str,
    reads_files: &Vec<String>,
    output_dir: &str,
    prefix_length: usize,
) {
    let singleton_kmers = load_kmer_db(bincode_file);
    let kmer_to_file = map_kmer_to_file(&singleton_kmers);
    let output_files = reads_files
        .par_iter()
        .map(|reads_file| classify_one(&singleton_kmers, &kmer_to_file, reads_file))
        .collect::<Vec<_>>();
    // Move output files to the output directory
    std::fs::create_dir_all(output_dir).expect("valid output directory");
    let new_output_files = output_files
        .iter()
        .map(|output_file| output_dir.to_string() + "/" + &output_file.split('/').last().unwrap())
        .collect::<Vec<_>>();
    for (output_file, new_output_file) in output_files.iter().zip(new_output_files.iter()) {
        std::fs::rename(output_file, new_output_file).expect("valid rename");
    }
    log::info!("Read classifications moved to `{}`", output_dir);

    // Collect the read classifications
    println!("Found {} read classifications", new_output_files.len());

    let dfs: Vec<Vec<ReadClassification>> = new_output_files
        .par_iter()
        .map(|rc| get_reads(rc, prefix_length))
        .collect();
    let mut all_reads = Vec::new();
    for df in dfs {
        all_reads.extend(df);
    }
    println!("Processed {} read classifications", all_reads.len());

    if !all_reads.is_empty() {
        let output_path = format!("{}.filtered.tsv", output_dir);
        let mut file = File::create(&Path::new(&output_path)).expect("Unable to create file");
        for read in all_reads {
            writeln!(file, "{}", read.join("\t")).expect("Unable to write row");
        }
    }
}

fn classify_one(
    singleton_kmers: &SingletonKmers,
    kmer_to_file: &HashMap<u64, usize>,
    reads_file: &str,
) -> String {
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

    output_file
}

fn get_reads(rc: &str, prefix_length: usize) -> Vec<ReadClassification> {
    let file = File::open(rc).expect("Unable to open file");
    let reader = BufReader::new(file);
    let mut filtered = Vec::new();

    for line in reader.lines() {
        let line = line.expect("Unable to read line");
        let row: Vec<String> = line.split('\t').map(String::from).collect();
        if row.len() < 4 {
            continue;
        }

        let kmers: i32 = row[2].parse().unwrap_or(0);
        let classification = &row[3];

        if classification.contains("Unclassified") {
            continue;
        }

        let parts: Vec<&str> = classification.splitn(2, ':').collect();
        let ab: Vec<&str> = parts[0].split(',').collect();

        if ab.len() < 2 || ab[0][..prefix_length] != ab[1][..prefix_length] {
            continue;
        }

        let (a, b) = if ab[0] < ab[1] {
            (ab[0], ab[1])
        } else {
            (ab[1], ab[0])
        };
        let scores: Vec<i32> = parts[1]
            .split(',')
            .map(|s| s.parse().unwrap_or(0))
            .collect();

        let mut new_row = row.clone();
        if kmers >= 300 && scores.iter().sum::<i32>() >= 50 && scores[1] >= 10 {
            new_row.push(format!("{}_{}", a, b));
            filtered.push(row);
        }
    }
    log::info!("Filtered {} reads", filtered.len());
    filtered
}
