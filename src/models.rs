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

impl ClassifyResults {
    pub fn tag(&self, fasta_files: &Vec<String>) -> String {
        let mut best_count = 0;
        let mut best_index = 0;
        let mut second_best_count = 0;
        let mut second_best_index = 0;
        let mut total = 0;
        for (i, &count) in self.scaled_counts.iter().enumerate() {
            total += count;
            if count > best_count {
                second_best_count = best_count;
                second_best_index = best_index;
                best_count = count;
                best_index = i;
            } else if count > second_best_count {
                second_best_count = count;
                second_best_index = i;
            }
        }
        let half = total / 2;
        if best_count + second_best_count > half {
            format!(
                "{},{}:{},{}",
                fasta_files[best_index],
                fasta_files[second_best_index],
                best_count * 100 / total,
                second_best_count * 100 / total
            )
        } else {
            format!(
                "Unclassified:{},{}",
                best_count * 100 / total,
                second_best_count * 100 / total
            )
        }
    }
}
