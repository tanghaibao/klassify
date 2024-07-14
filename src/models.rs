use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct SingletonKmers {
    pub kmer_size: u8,
    pub fasta_files: Vec<String>,
    pub kmers: Vec<Vec<u64>>,
}
