use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct SingletonKmers {
    pub kmer_size: u8,
    pub fasta_files: Vec<String>,
    pub kmers: Vec<Vec<u64>>,
    pub scaling_factors: Vec<f64>,
}

impl SingletonKmers {
    pub fn n(&self) -> usize {
        self.fasta_files.len()
    }
}
