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

pub struct ClassifyResults {
    pub id: String,
    pub seq_len: usize,
    pub counts: Vec<i32>,
    pub scaled_counts: Vec<i32>,
}

impl std::fmt::Display for ClassifyResults {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{} ({}; {} kmers) =>\n {}\n {}",
            self.id,
            self.seq_len,
            self.counts.iter().sum::<i32>(),
            self.counts
                .iter()
                .map(|count| count.to_string())
                .collect::<Vec<_>>()
                .join("\t"),
            self.scaled_counts
                .iter()
                .map(|count| count.to_string())
                .collect::<Vec<_>>()
                .join("\t")
        )
    }
}
