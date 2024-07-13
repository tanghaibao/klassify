use clap::Parser;

#[derive(Parser, Debug)]
pub struct BuildArgs {
    /// Input FASTA files
    pub fasta_files: Vec<String>,
}

pub fn build(fasta_files: &Vec<String>) {
    println!("{:?}", fasta_files);
}
