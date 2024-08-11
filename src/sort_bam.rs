use crate::models::MAX_DE;

use clap::Parser;
use log;
use num_cpus;
use rust_htslib::bam;
use rust_htslib::bam::{record::Aux, Read};

#[derive(Debug, Parser)]
pub struct SortBamArgs {
    /// Input BAM file
    pub input: String,
    /// Output BAM file
    #[clap(short, long)]
    pub output: String,
    /// Maximum divergence
    #[clap(short, long, default_value_t = MAX_DE)]
    pub max_de: f32,
}

/// Sort a BAM file by divergence
pub fn sort_bam(input_bam: &str, output_bam: &str, max_de: f32) {
    let n_threads = num_cpus::get();
    let n_threads_read = n_threads * 4 / 5;
    let n_threads_write = n_threads - n_threads_read;
    let mut bam = bam::Reader::from_path(input_bam).unwrap();
    bam.set_threads(n_threads).unwrap();
    log::info!(
        "Reading BAM file `{}` (n_threads_read={}, n_threads_write={})",
        input_bam,
        n_threads_read,
        n_threads_write
    );
    let header = bam::Header::from_template(bam.header());
    let mut out_bam = bam::Writer::from_path(output_bam, &header, bam::Format::Bam).unwrap();
    out_bam.set_threads(0).unwrap();
    let mut total_reads = 0;
    let mut retained_reads = 0;
    for r in bam.records() {
        total_reads += 1;
        if total_reads % 1_000_000 == 0 {
            log::info!(
                "Processed {} reads, retained {} reads ({}%)",
                total_reads,
                retained_reads,
                retained_reads * 100 / total_reads
            );
        }
        let rec = r.unwrap();
        let aux_cb = rec.aux(b"de").ok();
        if aux_cb.is_none() {
            continue;
        }
        let de = match aux_cb {
            Some(Aux::Float(de)) => de,
            _ => 0.0,
        };
        if de < max_de {
            out_bam.write(&rec).unwrap();
            retained_reads += 1;
        }
    }
    log::info!(
        "Sorted {} ({}% of {}) reads to `{}`",
        retained_reads,
        retained_reads * 100 / total_reads,
        total_reads,
        output_bam
    );
}
