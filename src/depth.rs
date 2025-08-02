use perbase_lib::par_granges::RegionProcessor;
use rust_htslib::bam::{self, Read};
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
        let mut rdr = bam::IndexedReader::from_path(&self.bam).unwrap();
        let hdr = rdr.header().to_owned();
        rdr.fetch((tid, start, stop)).unwrap();

        let mut bins = Vec::new();
        let mut sum = 0u64;
        let mut b_start = start;

        for p in rdr.pileup() {
            let pile = p.unwrap();
            let pos = pile.pos();
            // advance to the next bin if necessary
            while pos >= b_start + self.bin_size {
                bins.push(BedRecord {
                    chrom: String::from_utf8_lossy(hdr.tid2name(tid)).into(),
                    start: b_start,
                    end: b_start + self.bin_size,
                    depth: sum as f64 / self.bin_size as f64,
                });
                b_start += self.bin_size;
                sum = 0;
            }
            sum += pile.depth() as u64;
        }
        // flush the last (possibly partial) bin
        let end = (b_start + self.bin_size).min(stop);
        let count = end - b_start;
        if count > 0 {
            bins.push(BedRecord {
                chrom: String::from_utf8_lossy(hdr.tid2name(tid)).into(),
                start: b_start,
                end,
                depth: sum as f64 / count as f64,
            });
        }
        bins
    }
}
