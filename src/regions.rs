use crate::models::prefix;
use clap::Parser;
use log;
use perbase_lib::{
    par_granges::{self, RegionProcessor},
    position::pileup_position::PileupPosition,
    read_filter::ReadFilter,
};
use rayon::prelude::*;
use rust_htslib::bam::{self, pileup::Alignment, record::Record, Read};
use std::path::PathBuf;

// To use ParGranges you will need to implement a [`RegionProcessor`](par_granges::RegionProcessor),
// which requires a single method [`RegionProcessor::process_region`](par_granges::RegionProcessor::process_region)
// and an associated type P, which is the type of the values returned in the Vec by
// `process_region`. The returned `P` objects will be kept in order and accessible on the
// receiver channel returned by the `[ParGranges::process`](par_granges::ParGranges::process) method.
struct BasicProcessor<F: ReadFilter> {
    // An indexed bamfile to query for the region we were passed
    bamfile: PathBuf,
    // This is an object that implements `position::ReadFilter` and will be applied to
    // each read
    read_filter: F,
}

// A struct that holds the filter values that will be used to implement `ReadFilter`
struct BasicReadFilter {
    include_flags: u16,
    exclude_flags: u16,
    min_mapq: u8,
}

// The actual implementation of `ReadFilter`
impl ReadFilter for BasicReadFilter {
    // Filter reads based SAM flags and mapping quality, true means pass
    #[inline]
    fn filter_read(&self, read: &Record, _alignment: Option<&Alignment>) -> bool {
        let flags = read.flags();
        (!flags) & &self.include_flags == 0
            && flags & &self.exclude_flags == 0
            && &read.mapq() >= &self.min_mapq
    }
}

// Implementation of the `RegionProcessor` trait to process each region
impl<F: ReadFilter> RegionProcessor for BasicProcessor<F> {
    type P = PileupPosition;

    // This function receives an interval to examine.
    fn process_region(&self, tid: u32, start: u32, stop: u32) -> Vec<Self::P> {
        let mut reader = bam::IndexedReader::from_path(&self.bamfile).expect("Indexed reader");
        let header = reader.header().to_owned();
        // fetch the region
        reader.fetch((tid, start, stop)).expect("Fetched ROI");
        // Walk over pileups
        let result: Vec<PileupPosition> = reader
            .pileup()
            .flat_map(|p| {
                let pileup = p.expect("Extracted a pileup");
                // Verify that we are within the bounds of the chunk we are iterating on
                // Since pileup will pull reads that overhang edges.
                if pileup.pos() >= start && pileup.pos() < stop {
                    Some(PileupPosition::from_pileup(
                        pileup,
                        &header,
                        &self.read_filter,
                        None,
                    ))
                } else {
                    None
                }
            })
            .collect();
        result
    }
}

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
    // Create the read filter
    let read_filter = BasicReadFilter {
        include_flags: 0,
        exclude_flags: 3848,
        min_mapq: 20,
    };

    // Create the region processor
    let basic_processor = BasicProcessor {
        bamfile: bam_file.into(),
        read_filter: read_filter,
    };

    // Create a par_granges runner
    let par_granges_runner = par_granges::ParGranges::new(
        PathBuf::from(bam_file), // pass in bam
        None,                    // optional ref fasta
        None,                    // bedfile to narrow regions
        None,                    // optional bcf/vcf file to specify positions of interest
        false,                   // Merge any overlapping regions in the BED file
        None,                    // optional allowed number of threads, defaults to max
        None,                    // optional chunksize modification
        None, // optional modifier on the size of the channel for sending Positions
        basic_processor,
    );

    // Run the processor
    let receiver = par_granges_runner.process().unwrap();
    // let file_prefix = prefix(bam_file);
    // let output_bed = file_prefix + ".regions.bed";
    // Pull the in-order results from the receiver channel
    receiver.into_iter().for_each(|p: PileupPosition| {
        // Note that the returned values are required to be `serde::Serialize`, so more fancy things
        // than just debug printing are doable.
        println!("{:?}", p);
    });
}