use crate::models::{need_update, prefix_until_dot, sh};
use clap::Parser;
use csv::ReaderBuilder;
use flate2;
use log;
use rust_htslib::bam;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

/// Discrete bin size to contract regions
const BINSIZE: u32 = 10_000;
/// Chain distance to merge regions
const CHAIN_DISTANCE: u32 = 2 * BINSIZE;

#[derive(Parser, Debug)]
pub struct RegionsArgs {
    /// BAM files
    pub bam_files: Vec<String>,
    /// Retain only chimeras between chromosomes, must contain "Chr" and "chr"
    #[clap(short, long, default_value_t = false)]
    pub no_chr_only: bool,
}

#[derive(Debug)]
struct BedRecord {
    chrom: String,
    start: u32,
    end: u32,
    depth: String,
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
pub fn regions(bam_files: &Vec<String>, chr_only: bool) {
    let mut bed_files = Vec::new();
    for bam_file in bam_files {
        let bed_file = if bam_file.ends_with(".bam") {
            regions_one(bam_file)
        } else {
            bam_file.clone()
        };
        bed_files.push(bed_file);
    }

    // Perform the depth analysis
    process_bedfiles(bed_files, chr_only);
}

/// Prepare one BAM file and generate depths for each bin
fn regions_one(bam_file: &str) -> String {
    // Check if BAM index exists
    let bam_index = bam_file.to_string() + ".bai";
    if !std::path::Path::new(&bam_index).exists() {
        bam::index::build(bam_file, None, bam::index::Type::Bai, 1).unwrap();
        log::info!("Built index for `{}`", bam_file);
    }
    let mosdepth_out = prefix_until_dot(bam_file) + ".mosdepth";
    let mosdepth_cmd = format!(
        "mosdepth -t 8 -n --by {} {} {}",
        BINSIZE, mosdepth_out, bam_file
    );
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
    mosdepth_bed
}

/// Load BED file into a bunch of records
fn load_bed(bed: &str) -> Vec<BedRecord> {
    let file = BufReader::new(flate2::read::MultiGzDecoder::new(File::open(bed).unwrap()));
    let mut rdr = ReaderBuilder::new().delimiter(b'\t').from_reader(file);
    let mut records = Vec::new();

    for record in rdr.records() {
        let record = record.unwrap();
        records.push(BedRecord::from_csv_record(&record));
    }

    records
}

/// Process F1 and parent BED files to generate candidate regions.
fn process_bedfiles(bed_files: Vec<String>, chr_only: bool) -> HashMap<String, i32> {
    let child_bed = &bed_files[0];
    let parent1_bed = &bed_files[1];
    let parent2_bed = if bed_files.len() == 3 {
        Some(bed_files[2].clone())
    } else {
        None
    };

    let child_records = load_bed(child_bed);
    let parent1_records = load_bed(parent1_bed);
    let parent2_records = if let Some(bed) = parent2_bed {
        Some(load_bed(&bed))
    } else {
        None
    };

    let mut regions: BTreeMap<String, Vec<(u32, u32, f64)>> = BTreeMap::new();

    for (i, child_record) in child_records.iter().enumerate() {
        let child_depth = child_record.depth.parse::<f64>().unwrap();
        let parent1_depth = parent1_records[i].depth.parse::<f64>().unwrap();
        let parent2_depth = if let Some(ref parent2_recs) = parent2_records {
            parent2_recs[i].depth.parse::<f64>().unwrap()
        } else {
            0.0
        };

        let depth_ratio = child_depth / (parent1_depth + parent2_depth + 1.0);

        regions
            .entry(child_record.chrom.clone())
            .or_default()
            .push((child_record.start, child_record.end, depth_ratio));
    }

    let mut d = Vec::new();
    let mut selected = Vec::new();

    for (chrom, data) in &mut regions {
        if chr_only && !chrom.contains("Chr") && !chrom.contains("chr") {
            continue;
        }

        data.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap());

        let chrom_selected: Vec<_> = data
            .iter()
            .filter(|&&(_, _, depth)| depth >= 5.0 && depth <= 100.0)
            .map(|&(start, end, depth)| (chrom.clone(), start, end, format!("{}", depth.round())))
            .collect();

        selected.extend_from_slice(&chrom_selected);
        let regions_str = chrom_selected
            .iter()
            .map(|(chrom, start, end, depth)| {
                format!("{}:{}-{}:{}", chrom, start, end, depth.clone())
            })
            .collect::<Vec<_>>()
            .join(",");

        d.push((chrom.clone(), regions_str));
    }

    let prefix = Path::new(child_bed).file_stem().unwrap().to_str().unwrap();
    let poi_tsv = format!("{}.poi.tsv", prefix);

    let mut kf_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&poi_tsv)
        .unwrap();

    kf_writer.write_record(&["Chrom", "Regions"]).unwrap();

    for (chrom, regions_str) in d {
        kf_writer.write_record(&[chrom, regions_str]).unwrap();
    }

    log::info!("Points of interests written to `{}`", poi_tsv);

    // Merge regions that are close to each other
    selected.sort_by_key(|k| (k.0.clone(), k.1));

    let mut merged = Vec::new();

    for i in 0..selected.len() {
        if i == 0 {
            merged.push(selected[i].clone());
            continue;
        }

        let prev = merged.last_mut().unwrap();
        let cur = &selected[i];

        if prev.0 == cur.0 && prev.2 + CHAIN_DISTANCE >= cur.1 {
            prev.2 = prev.2.max(cur.2);
            prev.3 = format!("{},{}", prev.3, cur.3);
        } else {
            merged.push(cur.clone());
        }
    }

    // Write the merged regions to a file
    let regions_file = format!("{}.regions.tsv", prefix);
    let mut counter = HashMap::new();

    let mut regions_writer = BufWriter::new(File::create(&regions_file).unwrap());

    for (chrom, start, end, score) in &merged {
        writeln!(regions_writer, "{}:{}-{}\t{}", chrom, start, end, score).unwrap();
        *counter.entry(chrom[..2].to_string()).or_insert(0) += 1;
    }

    log::info!("Merged regions written to `{}`", regions_file);
    log::info!("Region counts: {:?}", counter);
    counter
}
