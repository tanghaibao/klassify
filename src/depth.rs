use perbase_lib::par_granges::RegionProcessor;
use serde::Serialize;
use std::path::PathBuf;

/// One BED-style line per fixed-size bin
#[derive(Debug, Serialize)]
pub(crate) struct BedRecord {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub depth: f64,
}

impl BedRecord {
    pub(crate) fn from_csv_record(record: &csv::StringRecord) -> BedRecord {
        BedRecord {
            chrom: record[0].to_string(),
            start: record[1].parse().unwrap(),
            end: record[2].parse().unwrap(),
            depth: record[3].parse().unwrap(),
        }
    }
}

/// Processor that collapses pile-ups into fixed-width windows
pub(crate) struct WindowProcessor {
    pub bam: PathBuf,
    pub bin_size: u32,
}
impl RegionProcessor for WindowProcessor {
    type P = BedRecord;

    fn process_region(&self, tid: u32, start: u32, stop: u32) -> Vec<Self::P> {
        use rust_htslib::bam::{self, Read};

        // BAM reader
        let mut rdr = bam::IndexedReader::from_path(&self.bam).unwrap();
        rdr.fetch((tid, start, stop)).unwrap();
        let chrom = String::from_utf8_lossy(rdr.header().tid2name(tid)).into_owned();

        // Bin bookkeeping
        let bin = self.bin_size;
        let n_bins = (stop - start).div_ceil(bin) as usize; // ceil div
        let mut sum = vec![0u64; n_bins]; // cumulative depth per bin
        let mut cov = vec![0u32; n_bins]; // # positions seen in bin

        // Stream pile-up once
        for pile in rdr.pileup() {
            let pile = pile.unwrap();
            let idx = ((pile.pos() - start) / bin) as usize;
            sum[idx] += pile.depth() as u64;
            cov[idx] += 1;
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
