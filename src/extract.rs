use crate::models::prefix;

use clap::Parser;
use csv::ReaderBuilder;
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

fn extract_one(read_map: &HashMap<String, String>, fasta_file: &str) {
    let mut reader = parse_fastx_file(fasta_file).expect("valid reads file");
    let file_prefix = prefix(fasta_file);
    let output_file = file_prefix.to_string() + ".extracted.fasta";
    let mut writer = std::fs::File::create(&output_file).unwrap();
    while let Some(record) = reader.next() {
        let record = record.expect("valid record");
        let seq = record.normalize(false);
        let id = String::from_utf8(record.id().to_vec())
            .unwrap()
            .split_whitespace()
            .next()
            .unwrap()
            .to_string();
        if let Some(label) = read_map.get(&id) {
            writeln!(
                writer,
                ">{}\n{}",
                label,
                String::from_utf8(seq.to_vec()).unwrap()
            )
            .unwrap();
        }
    }
}

pub fn extract(reads_tsv: &str, fasta_files: &Vec<String>) {
    let read_map = get_read_ids(reads_tsv);
    fasta_files.par_iter().for_each(|fasta_file| {
        extract_one(&read_map, fasta_file);
    });
}
