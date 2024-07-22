use crate::models::prefix;
use clap::Parser;
use log;
use perbase_lib::{
    par_granges::{self, RegionProcessor},
    position::{
        range_positions::{BedFormatRangePositions, RangePositions},
        Position,
    },
    read_filter::{DefaultReadFilter, ReadFilter},
    utils,
};
use rayon::prelude::*;
use rust_htslib::bam::{self, ext::BamRecordExtensions, Read};
use smartstring::alias::String;
use std::path::PathBuf;

/// Holds the info needed for [par_io::RegionProcessor] implementation
struct OnlyDepthProcessor<F: ReadFilter> {
    /// path to indexed BAM/CRAM
    reads: PathBuf,
    /// Indicate whether or not to keep positions with 0 depth
    keep_zeros: bool,
    /// implementation of [position::ReadFilter] that will be used
    read_filter: F,
}

impl<F: ReadFilter> OnlyDepthProcessor<F> {
    /// Sum the counts within the region to get the depths at each RangePosition
    #[inline]
    fn sum_counter(
        &self,
        counter: Vec<i32>,
        contig: &str,
        region_start: u32,
    ) -> Vec<RangePositions> {
        let mut sum: i32 = 0;
        let mut results = vec![];
        for (i, count) in counter.iter().enumerate() {
            sum += count;
            if sum != 0 || self.keep_zeros {
                let mut pos = RangePositions::new(String::from(contig), region_start + i as u32);
                pos.depth = u32::try_from(sum).expect("All depths are positive");
                pos.end = region_start + i as u32 + 1;
                results.push(pos);
            }
        }
        results
    }

    fn process_region_fast(&self, tid: u32, start: u32, stop: u32) -> Vec<RangePositions> {
        // Create a reader
        let mut reader =
            bam::IndexedReader::from_path(&self.reads).expect("Indexed Reader for region");

        let header = reader.header().to_owned();
        // fetch the region of interest
        reader.fetch((tid, start, stop)).expect("Fetched a region");

        let mut counter: Vec<i32> = vec![0; (stop - start) as usize];
        // Walk over each read, counting the starts and ends
        for record in reader
            .rc_records()
            .map(|r| r.expect("Read record"))
            .filter(|read| self.read_filter.filter_read(&read, None))
        {
            let rec_start = u32::try_from(record.reference_start()).expect("check overflow");
            let rec_stop = u32::try_from(record.reference_end()).expect("check overflow");

            // rectify start / stop with region boundaries
            // NB: impossible for rec_start > start since this is from fetch and we aren't splitting bam
            let adjusted_start = if rec_start < start {
                0
            } else {
                (rec_start - start) as usize
            };

            let mut dont_count_stop = false; // set this flag if this interval extends past the end of our region
            let adjusted_stop = if rec_stop >= stop {
                dont_count_stop = true;
                counter.len() - 1
            } else {
                (rec_stop - start) as usize
            };

            counter[adjusted_start] += 1;
            // Check if end of interval extended past region end
            if !dont_count_stop {
                counter[adjusted_stop] -= 1;
            }
        }

        // Sum the counter and merge same-depth ranges of positions
        let contig = std::str::from_utf8(header.tid2name(tid)).unwrap();
        self.sum_counter(counter, contig, start)
    }
}

/// Implement [par_io::RegionProcessor] for [SimpleProcessor]
impl<F: ReadFilter> RegionProcessor for OnlyDepthProcessor<F> {
    /// Objects of [position::Position] will be returned by each call to [SimpleProcessor::process_region]
    type P = RangePositions;

    /// Process a region by fetching it from a BAM/CRAM, getting a pileup, and then
    /// walking the pileup (checking bounds) to create Position objects according to
    /// the defined filters
    fn process_region(&self, tid: u32, start: u32, stop: u32) -> Vec<RangePositions> {
        self.process_region_fast(tid, start, stop)
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
    let read_filter = DefaultReadFilter::new(0, 512, 0);

    // Create the region processor
    let depth_processor = OnlyDepthProcessor {
        reads: PathBuf::from(bam_file),
        keep_zeros: true,
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
        depth_processor,
    );

    // Run the processor
    let receiver = par_granges_runner.process().unwrap();
    let file_prefix = prefix(bam_file);
    let output_bed = file_prefix + ".regions.bed.gz";
    let mut writer = utils::get_writer(&Some(output_bed), true, false, 4, 2).unwrap();
    // Pull the in-order results from the receiver channel
    for pos in receiver.into_iter() {
        writer
            .serialize(BedFormatRangePositions::from(pos))
            .unwrap();
    }
    writer.flush().unwrap();
}
