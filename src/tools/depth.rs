use perbase_lib::par_granges::{ParGranges, RegionProcessor};
use perbase_lib::position::pileup_position::PileupPosition;
use perbase_lib::read_filter::ReadFilter;
use rust_htslib::bam::pileup::Alignment;
use rust_htslib::bam::Record;
use rust_htslib::bam::{self, Read};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

/// A runner that processes BAM files with a basic read filter
pub type Runner = ParGranges<WindowProcessor<BasicReadFilter>>;

/// One BED-style line per fixed-size bin
#[derive(Debug, Deserialize, Serialize)]
pub struct BedRecord {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub depth: f64,
}

/// A struct that holds the filter values that will be used to implement `ReadFilter`
pub struct BasicReadFilter {
    include_flags: u16,
    exclude_flags: u16,
    min_mapq: u8,
}

/// The actual implementation of `ReadFilter`
impl ReadFilter for BasicReadFilter {
    // Filter reads based SAM flags and mapping quality, true means pass
    #[inline]
    fn filter_read(&self, read: &Record, _: Option<&Alignment>) -> bool {
        let flags = read.flags();
        (!flags) & self.include_flags == 0
            && flags & self.exclude_flags == 0
            && read.mapq() >= self.min_mapq
    }
}

/// Processor that collapses pile-ups into fixed-width windows
pub struct WindowProcessor<F: ReadFilter> {
    /// An indexed bamfile to query for the region we were passed
    pub bam: PathBuf,
    /// Size of the bins to collapse the region into
    pub bin_size: u32,
    /// Read filter to apply to the reads
    read_filter: F,
}

impl<F: ReadFilter> RegionProcessor for WindowProcessor<F> {
    type P = BedRecord;

    fn process_region(&self, tid: u32, start: u32, stop: u32) -> Vec<Self::P> {
        // BAM reader
        let mut rdr = bam::IndexedReader::from_path(&self.bam).unwrap();
        rdr.fetch((tid, start, stop)).unwrap();
        let header = rdr.header().to_owned();
        let chrom = String::from_utf8_lossy(header.tid2name(tid)).into_owned();

        // Bin bookkeeping
        let bin = self.bin_size;
        let n_bins = (stop - start).div_ceil(bin) as usize; // ceil div
        let mut sum = vec![0u64; n_bins]; // cumulative depth per bin

        // Stream pile-up once
        for pile in rdr.pileup() {
            let pile = pile.unwrap();
            let pileup = PileupPosition::from_pileup(pile, &header, &self.read_filter, None);
            let pos = pileup.pos;
            if pos < start || pos >= stop {
                continue; // skip positions outside the region
            }
            let idx = ((pos - start) / bin) as usize;
            sum[idx] += pileup.depth as u64; // accumulate depth
        }

        // Convert to BED records
        (0..n_bins)
            .map(|i| {
                let b_start = start + (i as u32) * bin;
                let b_end = (b_start + bin).min(stop); // last bin may be short
                let width = (b_end - b_start) as f64; // denom
                BedRecord {
                    chrom: chrom.clone(),
                    start: b_start,
                    end: b_end,
                    depth: sum[i] as f64 / width, // mean depth
                }
            })
            .collect()
    }
}

/// Create a new runner for processing BAM files with a basic read filter
pub fn get_runner(bam: &str, bin_size: u32) -> Runner {
    let bam = PathBuf::from(bam);
    // Create the read filter
    let read_filter = BasicReadFilter {
        include_flags: 0,
        exclude_flags: 1796, // mosdepth default: unmapped, secondary, QC failed, duplicate
        min_mapq: 0,
    };

    // Create the region processor
    let proc = WindowProcessor {
        bam: bam.clone(),
        bin_size,
        read_filter,
    };
    ParGranges::new(
        bam, None, None, None,  // no fasta, no VCF/BED constraints
        false, // merge intervals
        None, None, None, // default threads/chunking
        proc,
    )
}
