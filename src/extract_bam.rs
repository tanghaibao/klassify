use clap::Parser;
use log;
use rust_htslib::bam::{IndexedReader, Read};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};

#[derive(Parser, Debug)]
pub struct ExtractBamArgs {
    /// Regions that contain the reads
    /// Format: "chr1:1-100", one region per line
    pub regions_file: String,
    /// BAM file to extract reads
    pub bam_file: String,
}

/// Extract reads from a BAM file
pub fn extract_bam(region_file: &str, bam_file: &str) {
    log::info!("Extracting reads from BAM file");
    let regions = std::fs::read_to_string(region_file).expect("valid region file");
    let regions: Vec<&str> = regions.lines().collect();
    let mut bam = IndexedReader::from_path(bam_file).expect("valid BAM file");
    let output_file = region_file.to_string() + ".extracted.fasta";
    let mut writer = BufWriter::new(File::create(&output_file).unwrap());
    let mut seen = HashSet::new();
    for region in regions {
        let region = region.split(':').collect::<Vec<&str>>();
        let chrom = region[0];
        let region = region[1].split('-').collect::<Vec<&str>>();
        let start = region[0].parse::<i32>().expect("valid start position");
        let end = region[1].parse::<i32>().expect("valid end position");
        let _ = bam.fetch((chrom, start, end)).expect("valid region");
        for record in bam.records() {
            let record = record.expect("valid record");
            if record.is_unmapped()
                || record.is_secondary()
                || record.is_quality_check_failed()
                || record.is_duplicate()
                || record.is_supplementary()
            {
                continue;
            }
            let id = String::from_utf8(record.qname().to_vec()).expect("valid read ID");
            if seen.contains(&id) {
                log::info!("Skipping duplicate read {}", id);
                continue;
            }
            seen.insert(id.clone());
            let seq = record.seq();
            writeln!(
                writer,
                ">{}\n{}",
                id,
                String::from_utf8(seq.as_bytes()).unwrap()
            )
            .expect("write to FASTA");
        }
    }
}
