use crate::models::SingletonKmers;
use bincode::deserialize_from;
use clap::Parser;
use log;
use std::fs::File;

#[derive(Parser, Debug)]
pub struct ClassifyArgs {
    /// Bincode file
    pub bincode_file: String,
    /// Read file to classify
    pub reads_file: String,
}

pub fn classify(bincode_file: &str, reads_file: &str) {
    let _: SingletonKmers = deserialize_from(File::open(bincode_file).unwrap()).unwrap();
    log::info!("Loaded singleton kmers");
}
