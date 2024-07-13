use clap::Parser;
use needletail::{parse_fastx_file, Sequence};
use std::collections::HashSet;

const KMER_SIZE: u8 = 24;

#[derive(Parser, Debug)]
pub struct BuildArgs {
    /// Input FASTA files
    pub fasta_files: Vec<String>,
}

pub fn build(fasta_files: &Vec<String>) {
    for fasta_file in fasta_files {
        get_kmers(fasta_file);
    }
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
    log::info!("{} kmers found", kmer_set.len());
    kmer_set
}
