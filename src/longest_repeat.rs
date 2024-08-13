use crate::models::DEFAULT_MIN_REPEAT_LENGTH;

use bio::data_structures::suffix_array::{lcp, suffix_array};
use clap::Parser;
use log;
use needletail::{parse_fastx_file, Sequence};

#[derive(Parser, Debug)]
pub struct LongestRepeatArgs {
    /// FASTA file to extract reads
    pub fasta_file: String,
    /// Minimum length of the repeated substring
    #[clap(short, long, default_value_t = DEFAULT_MIN_REPEAT_LENGTH)]
    pub min_length: isize,
}

/// Find the longest repeated substring in a set of sequences
pub fn longest_repeat(fasta_file: &str, min_length: isize) {
    let mut reader = parse_fastx_file(fasta_file).expect("valid reads file");
    let mut text = Vec::new();
    let mut seq_names = Vec::new();
    let mut seq_starts = Vec::new();
    while let Some(record) = reader.next() {
        let record = record.expect("valid record");
        let seq = record.normalize(false);
        seq_names.push(String::from_utf8(record.id().to_vec()).unwrap());
        seq_starts.push(text.len());
        text.extend(seq.into_owned());
        text.push(b'$');
    }
    log::info!("Read {} sequences from FASTA file", seq_starts.len());

    let sa = suffix_array(&text);
    log::info!("Computed suffix array: {}", sa.len());
    let lcp = lcp(&text, &sa);
    log::info!("Computed LCP array: {}", lcp.len());
    for (index, lc) in lcp.iter().enumerate() {
        if lc >= min_length {
            // let seq_index = seq_starts.binary_search(&sa[index]).unwrap();
            let text_index = sa[index];
            let next_text_index = sa[index + 1];
            println!(
                "lcp={} => {}… {}…",
                lc,
                String::from_utf8(text[text_index..text_index + 20].to_vec()).unwrap(),
                String::from_utf8(text[next_text_index..next_text_index + 20].to_vec()).unwrap()
            );
        }
    }
}
