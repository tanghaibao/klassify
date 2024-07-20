use crate::models::prefix;

use clap::Parser;
use csv::ReaderBuilder;
use log;
use needletail::{parse_fastx_file, Sequence};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write;

#[derive(Parser, Debug)]
pub struct ExtractArgs {
    /// Filtered reads TSV file
    pub reads_tsv: String,
    /// FASTA files to extract kmers
    pub fasta_files: Vec<String>,
    /// Output file
    #[clap(short, long, default_value = "extracted.fasta")]
    pub output_file: String,
}

fn get_read_ids(reads_tsv: &str) -> HashMap<String, String> {
    let mut read_map = HashMap::new();
    let reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(reads_tsv)
        .expect("valid TSV file");
    for result in reader.into_records() {
        let record = result.expect("valid record");
        let read_id = record.get(0).expect("valid read ID").to_string();
        let label = record
            .get(record.len() - 1)
            .expect("valid label")
            .to_string();
        let new_read_id = label + "_" + &read_id;
        read_map.insert(read_id, new_read_id);
    }
    read_map
}

fn extract_one(read_map: &HashMap<String, String>, fasta_file: &str) -> String {
    let mut reader = parse_fastx_file(fasta_file).expect("valid reads file");
    let file_prefix = prefix(fasta_file);
    let output_file = file_prefix.to_string() + ".extracted.fasta";
    let mut writer = std::fs::File::create(&output_file).unwrap();
    let mut total = 0;
    while let Some(record) = reader.next() {
        let record = record.expect("valid record");
        let seq = record.normalize(false);
        let id = String::from_utf8(record.id().to_vec())
            .unwrap()
            .split_whitespace()
            .next()
            .unwrap()
            .to_string();
        if let Some(new_read_id) = read_map.get(&id) {
            writeln!(
                writer,
                ">{}\n{}",
                new_read_id,
                String::from_utf8(seq.to_vec()).unwrap()
            )
            .unwrap();
            total += 1;
        }
    }
    log::info!("Extracted {} reads from `{}`", total, fasta_file);
    output_file
}

/// Extract reads from FASTA/FASTQ files
pub fn extract(reads_tsv: &str, fasta_files: &Vec<String>, output_file: &str) {
    let read_map = get_read_ids(reads_tsv);
    let output_files = fasta_files
        .par_iter()
        .map(|fasta_file| extract_one(&read_map, fasta_file))
        .collect::<Vec<_>>();
    // Merge the output files
    let mut writer = std::fs::File::create(output_file).unwrap();
    for output_file in output_files.iter() {
        let mut reader = std::fs::File::open(&output_file).unwrap();
        std::io::copy(&mut reader, &mut writer).expect("valid copy");
    }
    // Cleanup the output files
    for output_file in output_files {
        std::fs::remove_file(output_file).expect("valid remove");
    }
    log::info!("Extracted reads written to `{}`", output_file);
}
