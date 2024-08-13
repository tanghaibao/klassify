use crate::models::DEFAULT_MIN_REPEAT_LENGTH;

use bio::data_structures::suffix_array::{lcp, suffix_array};
use clap::Parser;
use log;
use needletail::{parse_fastx_file, Sequence};
use std::fmt::Display;

const TERMINATOR: u8 = b'$';

#[derive(Parser, Debug)]
pub struct LongestRepeatArgs {
    /// FASTA file to extract reads
    pub fasta_file: String,
    /// Minimum length of the repeated substring
    #[clap(short, long, default_value_t = DEFAULT_MIN_REPEAT_LENGTH)]
    pub min_length: isize,
}

/// A range of positions in a sequence
struct SeqRange {
    name: String,
    start: usize,
    end: usize,
}

impl Display for SeqRange {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}-{}", self.name, self.start, self.end)
    }
}

/// Container for multiple sequences, concatenated together
struct MultiSequence {
    text: Vec<u8>,
    names: Vec<String>,
    starts: Vec<usize>,
}

impl MultiSequence {
    /// Get the sequence and its start and end positions in the concatenated text
    fn slice(&self, start: usize, end: usize) -> Option<SeqRange> {
        let start_index = self.bisect_left(start);
        let end_index = self.bisect_left(end);
        if start_index != end_index {
            return None;
        }
        let name = self.names[start_index].clone();
        let seq_start = self.starts[start_index];
        let new_start = start - seq_start;
        let mut new_end = end - seq_start;
        if self.text[new_end] == TERMINATOR {
            new_end -= 1;
        }
        Some(SeqRange {
            name,
            start: new_start,
            end: new_end,
        })
    }

    /// Find index of the sequence containing the given position
    #[inline]
    fn bisect_left(&self, pos: usize) -> usize {
        self.starts.binary_search(&pos).unwrap_or_else(|x| x - 1)
    }
}

/// Find the longest repeated substring in a set of sequences
pub fn longest_repeat(fasta_file: &str, min_length: isize) {
    let mut reader = parse_fastx_file(fasta_file).expect("valid reads file");
    let mut text = Vec::new();
    let mut names = Vec::new();
    let mut starts = Vec::new();
    while let Some(record) = reader.next() {
        let record = record.expect("valid record");
        let seq = record.normalize(false);
        let name = String::from_utf8(record.id().to_vec()).unwrap();
        names.push(name);
        starts.push(text.len());
        text.extend(seq.into_owned());
        text.push(TERMINATOR);
    }
    let ms = MultiSequence {
        text,
        names,
        starts,
    };
    log::info!("Read {} sequences from FASTA file", ms.names.len());

    let sa = suffix_array(&ms.text);
    log::info!("Computed suffix array: {}", sa.len());
    let lcp = lcp(&ms.text, &sa);
    log::info!("Computed LCP array: {}", lcp.len());
    const MAX_PRINT: usize = 16;
    for (index, lc) in lcp.iter().enumerate() {
        if lc < min_length {
            continue;
        }
        let start = sa[index - 1];
        let end = start + lc as usize;
        let next_start = sa[index];
        let next_end = next_start + lc as usize;
        let seqrange = ms.slice(start, end);
        let next_seqrange = ms.slice(next_start, next_end);
        if seqrange.is_none() || next_seqrange.is_none() {
            continue;
        }
        let seqrange = seqrange.unwrap();
        let next_seqrange = next_seqrange.unwrap();
        println!("{} => {}", seqrange, next_seqrange);
        let t1: String = String::from_utf8(ms.text[start..end].to_vec()).unwrap();
        let t2 = String::from_utf8(ms.text[next_start..next_end].to_vec()).unwrap();
        assert_eq!(t1, t2);
        println!(
            "lcp={} => {} … {} {} … {}",
            lc,
            t1[..MAX_PRINT].to_string(),
            t1[t1.len() - MAX_PRINT..].to_string(),
            t2[..MAX_PRINT].to_string(),
            t2[t2.len() - MAX_PRINT..].to_string()
        );
    }
}
