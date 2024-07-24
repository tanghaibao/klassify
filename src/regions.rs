use crate::models::{need_update, prefix_until_dot, sh};
use clap::Parser;
use csv::ReaderBuilder;
use log;
use rayon::prelude::*;
use rust_htslib::bam;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

#[derive(Parser, Debug)]
pub struct RegionsArgs {
    /// BAM files
    pub bam_files: Vec<String>,
}

#[derive(Debug)]
struct BedRecord {
    chrom: String,
    start: u32,
    end: u32,
    depth: f64,
}

impl BedRecord {
    fn from_csv_record(record: &csv::StringRecord) -> BedRecord {
        BedRecord {
            chrom: record[0].to_string(),
            start: record[1].parse().unwrap(),
            end: record[2].parse().unwrap(),
            depth: record[3].parse().unwrap(),
        }
    }
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
    let mosdepth_bed = mosdepth_out.to_string() + ".regions.bed.gz";
    if need_update(
        vec![bam_file.to_string()],
        vec![mosdepth_bed.to_string()],
        true,
    ) {
        let status = sh(&mosdepth_cmd);
        if status {
            log::info!("Generated depths for `{}`", bam_file);
        } else {
            log::error!("Failed to generate depths for `{}`", bam_file);
        }
    }
}
