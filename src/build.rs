use clap::Parser;

#[derive(Parser, Debug)]
pub struct BuildArgs {
    /// Input FASTA files
    #[clap(short, long)]
    pub input: Vec<String>,
}

pub fn build(args: &Vec<String>) {
    println!("{:?}", args);
}
