use crate::models::{prefix_until_dot, sh};
use clap::Parser;
use log;
use rayon::prelude::*;
use rust_htslib::bam;

#[derive(Parser, Debug)]
pub struct RegionsArgs {
    /// BAM files
    pub bam_files: Vec<String>,
}

/// Prepare BAM files and generate depths for each bin
pub fn regions(bam_files: &Vec<String>) {
    bam_files
        .par_iter()
        .for_each(|bam_file| regions_one(bam_file));
}

/// Prepare one BAM file and generate depths for each bin
fn regions_one(bam_file: &str) {
    // Check if BAM index exists
    let bam_index = bam_file.to_string() + ".bai";
    if !std::path::Path::new(&bam_index).exists() {
        bam::index::build(bam_file, None, bam::index::Type::Bai, 1).unwrap();
        log::info!("Built index for `{}`", bam_file);
    }
    let mosdepth_out = prefix_until_dot(bam_file) + ".mosdepth";
    let mosdepth_cmd = format!("mosdepth -t 8 -n --by 10000 {} {}", mosdepth_out, bam_file);
    let status = sh(&mosdepth_cmd);
    if status {
        log::info!("Generated depths for `{}`", bam_file);
    } else {
        log::error!("Failed to generate depths for `{}`", bam_file);
    }
}
