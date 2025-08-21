use crate::tools::graph_matching::greedy_max_weight_matching;
use clap::Parser;
use log::{info, warn};
use rust_htslib::bam::{self, record::Cigar::*, Read};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Default minimum read support for clustering
const DEFAULT_MIN_READ_SUPPORT: usize = 3;
/// Separator used in read IDs
const SEQID_MISMATCH: &str = "SEQID_MISMATCH";
/// Secondary match (not the longest span) in a BAM record
const SECONDARY_MATCH: &str = "SECONDARY_MATCH";

#[derive(Parser, Debug)]
pub struct ClusterPairsArgs {
    /// Input BAM containing alignments of split subreads ("orig|seqid|start-end")
    pub bam: String,
    /// Minimum read support required for regions and pairs
    #[clap(long, default_value_t = DEFAULT_MIN_READ_SUPPORT)]
    pub min_read_support: usize,
}

/// Maps from reads to regions, subreads, and matches.
pub struct ReadMaps {
    pub read_to_regions: HashMap<String, Vec<String>>,
    pub read_to_subreads: HashMap<String, HashSet<(String, char)>>,
    pub read_to_match: HashMap<String, String>,
}

/// Cluster pairs of subreads based on their alignment in a BAM file.
pub fn cluster_pairs(args: &ClusterPairsArgs) {
    // ------- Pass 1: track the longest reference span per subread (QNAME) -------
    let mut seen_max: HashMap<String, i64> = HashMap::new();
    {
        let mut r = bam::Reader::from_path(&args.bam).expect("open BAM");
        for rec in r.records() {
            let rec = rec.expect("record");
            if rec.is_unmapped() {
                continue;
            }
            let span = ref_span(&rec);
            let qn = qname(&rec);
            let e = seen_max.entry(qn).or_insert(0);
            if span > *e {
                *e = span;
            }
        }
    }

    // ------- Pass 2: filter primary (by max span) & seqid match -------------
    let mut counters: HashMap<&'static str, usize> = HashMap::new();
    let mut hits: Vec<Hit> = Vec::new();

    let mut rdr = bam::Reader::from_path(&args.bam).expect("open BAM");
    let hv = rdr.header().to_owned();

    for rec in rdr.records() {
        let rec = rec.expect("record");
        if rec.is_unmapped() {
            continue;
        }
        let qn = qname(&rec);
        let parts: Vec<&str> = qn.split('|').collect();
        if parts.len() < 3 {
            warn!(
                "QNAME `{}` not in `orig|seqid|start-end` format; dropped",
                qn
            );
            *counters.entry(SEQID_MISMATCH).or_insert(0) += 1;
            continue;
        }
        let expected_seqid = parts[1];

        let tid = rec.tid();
        if tid < 0 {
            continue;
        }
        let rname = String::from_utf8(hv.tid2name(tid as u32).to_vec()).unwrap();
        if rname != expected_seqid {
            *counters.entry(SEQID_MISMATCH).or_insert(0) += 1;
            continue;
        }

        let span = ref_span(&rec);
        if let Some(&mx) = seen_max.get(&qn) {
            if span != mx {
                *counters.entry(SECONDARY_MATCH).or_insert(0) += 1;
                continue;
            }
        }

        let start = rec.pos(); // 0-based
        let end = start + span;
        let strand = if rec.is_reverse() { '-' } else { '+' };

        hits.push(Hit {
            read_accn: qn.to_string(),
            read_root: parts[0].to_string(),
            ref_name: rname,
            start,
            end,
            strand,
        });
    }

    info!("Counters: {:?}", counters);
    info!("Total filtered: {}", hits.len());

    // ------- Cluster overlapping hits on the same reference -----------------
    hits.sort_by(|a, b| {
        (a.ref_name.as_str(), a.start, a.end).cmp(&(b.ref_name.as_str(), b.start, b.end))
    });

    let clusters = cluster_overlaps(hits);

    // ------- Build per-read maps from supported clusters --------------------
    let ReadMaps {
        read_to_regions,
        read_to_subreads,
        read_to_match,
    } = build_read_maps(clusters, args.min_read_support);

    // ------- Pair regions by reads with exactly two regions -----------------
    let pair_to_reads = build_pair_to_reads(&read_to_regions);

    // ------- Greedy max-weight matching (disjoint pairs) --------------------
    let selected_pairs = greedy_max_weight_matching(&pair_to_reads, args.min_read_support);

    // ------- Orient pairs using majority subread order ----------------------
    let filtered_pair_to_reads =
        orient_and_collect(selected_pairs, &pair_to_reads, &read_to_subreads);

    // ------- Emit table to stdout + paired regions file ---------------------
    print_table(&filtered_pair_to_reads, &read_to_match);

    let paired_regions_file = format!(
        "{}.paired.regions",
        Path::new(&args.bam)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("out")
    );
    let mut fw = BufWriter::new(File::create(&paired_regions_file).expect("create paired.regions"));
    for (ra, rb) in filtered_pair_to_reads.keys() {
        writeln!(fw, "{}", ra).unwrap();
        writeln!(fw, "{}", rb).unwrap();
    }
    info!("Paired regions written to `{}`", paired_regions_file);
}

// =========================
// Internal helpers (public for tests)
// =========================

/// Represents a hit in the BAM file, containing information about the read and its alignment.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Hit {
    pub read_accn: String, // "orig|seqid|start-end"
    pub read_root: String, // "orig"
    pub ref_name: String,  // reference (RNAME)
    pub start: i64,
    pub end: i64,
    pub strand: char, // '+'|'-'
}

/// Get the QNAME from a BAM record as a String.
fn qname(rec: &bam::Record) -> String {
    String::from_utf8(rec.qname().to_vec()).unwrap()
}

/// Calculate the reference span of a BAM record based on its CIGAR operations.
fn ref_span(rec: &bam::Record) -> i64 {
    rec.cigar()
        .iter()
        .map(|c| match *c {
            Match(l) | Del(l) | RefSkip(l) | Equal(l) | Diff(l) => l as i64,
            _ => 0,
        })
        .sum()
}

/// Calculate the median of a vector of i64 values.
fn median_i64(mut vals: Vec<i64>) -> i64 {
    if vals.is_empty() {
        return 0;
    }
    let mid = vals.len() / 2;
    vals.select_nth_unstable(mid);
    vals[mid]
}

/// Parse a read accession in the format "orig|seqid|start-end" into its components.
fn parse_accn(accn: &str) -> (&str, &str, &str) {
    let mut sp = accn.split('|');
    let read = sp.next().unwrap_or("");
    let seq = sp.next().unwrap_or("");
    let rng = sp.next().unwrap_or("");
    (read, seq, rng)
}

/// Cluster hits based on overlapping regions in the same reference.
pub fn cluster_overlaps(hits: Vec<Hit>) -> Vec<Vec<Hit>> {
    let mut clusters: Vec<Vec<Hit>> = Vec::new();
    for h in hits.into_iter() {
        if let Some(last) = clusters.last_mut() {
            if !last.is_empty()
                && last.last().unwrap().ref_name == h.ref_name
                && last.last().unwrap().end >= h.start
            {
                last.push(h);
                continue;
            }
        }
        clusters.push(vec![h]);
    }
    clusters
}

/// Build maps from reads to regions, subreads, and matches based on clusters of hits.
pub fn build_read_maps(clusters: Vec<Vec<Hit>>, min_read_support: usize) -> ReadMaps {
    let mut read_to_regions: HashMap<String, Vec<String>> = HashMap::new();
    let mut read_to_subreads: HashMap<String, HashSet<(String, char)>> = HashMap::new();
    let mut read_to_match: HashMap<String, String> = HashMap::new();

    for mb in clusters.into_iter() {
        if mb.len() < min_read_support {
            continue;
        }
        let mb_start = median_i64(mb.iter().map(|x| x.start).collect());
        let mb_end = median_i64(mb.iter().map(|x| x.end).collect());
        let region_name = format!("{}:{}-{}", mb[0].ref_name, mb_start, mb_end);

        for b in mb.into_iter() {
            let read_name = b.read_root.clone();
            read_to_regions
                .entry(read_name.clone())
                .or_default()
                .push(region_name.clone());
            read_to_subreads
                .entry(read_name.clone())
                .or_default()
                .insert((b.read_accn.clone(), b.strand));
            read_to_match.insert(
                b.read_accn.clone(),
                format!("{}:{}-{}:{}", b.ref_name, b.start, b.end, b.strand),
            );
        }
    }
    ReadMaps {
        read_to_regions,
        read_to_subreads,
        read_to_match,
    }
}

/// Build a mapping from pairs of regions to reads that support them.
pub fn build_pair_to_reads(
    read_to_regions: &HashMap<String, Vec<String>>,
) -> HashMap<(String, String), Vec<String>> {
    let mut pair_to_reads: HashMap<(String, String), Vec<String>> = HashMap::new();
    for (read, regions) in read_to_regions.iter() {
        if regions.len() == 2 {
            let (a, b) = if regions[0] <= regions[1] {
                (regions[0].clone(), regions[1].clone())
            } else {
                (regions[1].clone(), regions[0].clone())
            };
            pair_to_reads.entry((a, b)).or_default().push(read.clone());
        }
    }
    pair_to_reads
}
/// Orient pairs of reads based on their subread orientations and collect them.
pub fn orient_and_collect(
    selected_pairs: Vec<(String, String)>,
    pair_to_reads: &HashMap<(String, String), Vec<String>>,
    read_to_subreads: &HashMap<String, HashSet<(String, char)>>,
) -> BTreeMap<(String, String), Vec<(String, String)>> {
    let mut out: BTreeMap<(String, String), Vec<(String, String)>> = BTreeMap::new();

    for (ra, rb) in selected_pairs.into_iter() {
        let reads = &pair_to_reads[&(ra.clone(), rb.clone())];

        // (fa_seqid, fb_seqid) -> count
        let mut counter: HashMap<(String, String), usize> = HashMap::new();
        let mut gathered: Vec<(String, String)> = Vec::new();

        for read in reads.iter() {
            let subs = match read_to_subreads.get(read) {
                Some(s) if s.len() == 2 => s,
                _ => continue,
            };
            let mut it = subs.iter();
            let (mut fa, _sa) = it.next().unwrap().clone();
            let (mut fb, _sb) = it.next().unwrap().clone();

            let (_ra, _sa, _) = parse_accn(&fa);
            let (_rb, _sb2, fb_range) = parse_accn(&fb);
            if fb_range.starts_with("0-") {
                std::mem::swap(&mut fa, &mut fb);
            }
            let (_r1, fa_seqid, fa_rng) = parse_accn(&fa);
            let (_r2, fb_seqid, _fb_rng) = parse_accn(&fb);
            assert!(
                !fa_rng.ends_with("0-"),
                "left piece should not end with '0-'"
            );

            *counter
                .entry((fa_seqid.to_string(), fb_seqid.to_string()))
                .or_insert(0) += 1;
            gathered.push((fa, fb));
        }

        if counter.is_empty() {
            continue;
        }

        // Most common subread orientation vote (ra_reo, rb_reo)
        let ((ra_reo, rb_reo), _cnt) = counter.into_iter().max_by(|a, b| a.1.cmp(&b.1)).unwrap();

        // Ensure pair names match orientation (left starts with ra_reo)
        let mut oa = ra.clone();
        let mut ob = rb.clone();
        if !oa.starts_with(&ra_reo) {
            std::mem::swap(&mut oa, &mut ob);
        }

        // If still inconsistent, skip this pair (don't panic)
        if !(oa.starts_with(&ra_reo) && ob.starts_with(&rb_reo)) {
            warn!(
                "Orientation mismatch: pair=({},{}) vote=({},{}) — skipping",
                ra, rb, ra_reo, rb_reo
            );
            continue;
        }

        out.insert((oa, ob), gathered);
    }

    out
}

/// Print the filtered pairs and their reads in a tabular format.
fn print_table(
    filtered_pair_to_reads: &BTreeMap<(String, String), Vec<(String, String)>>,
    read_to_match: &HashMap<String, String>,
) {
    let header = [
        "Crossover ID",
        "Left",
        "Right",
        "Read Count",
        "Read ID",
        "Read Left",
        "Read Left Match",
        "Read Right",
        "Read Right Match",
    ];
    println!("{}", header.join("\t"));

    let mut cid = 0usize;
    for ((ra, rb), reads) in filtered_pair_to_reads.iter() {
        cid += 1;
        for (i, (fa, fb)) in reads.iter().enumerate() {
            let mut row: Vec<String> = Vec::new();
            if i == 0 {
                row.push(cid.to_string());
                row.push(ra.clone());
                row.push(rb.clone());
                row.push(reads.len().to_string());
            } else {
                row.extend([cid.to_string(), String::new(), String::new(), String::new()]);
            }
            let read_id = fa.split('|').next().unwrap_or("").to_string();
            let left_slice = fa.split('|').nth(2).unwrap_or("").to_string();
            let right_slice = fb.split('|').nth(2).unwrap_or("").to_string();
            let lm = read_to_match.get(fa).cloned().unwrap_or_default();
            let rm = read_to_match.get(fb).cloned().unwrap_or_default();

            row.push(read_id);
            row.push(left_slice);
            row.push(lm);
            row.push(right_slice);
            row.push(rm);
            println!("{}", row.join("\t"));
        }
    }
}

// =========================
// Tests
// =========================

#[cfg(test)]
mod tests {
    use super::*;

    fn mk_hit(root: &str, seq: &str, st: i64, en: i64, strand: char) -> Hit {
        Hit {
            read_accn: format!("{}|{}|{}-{}", root, seq, st, en),
            read_root: root.to_string(),
            ref_name: seq.to_string(),
            start: st,
            end: en,
            strand,
        }
    }

    #[test]
    fn test_cluster_overlaps_basic() {
        let hits = vec![
            mk_hit("r1", "A", 10, 20, '+'),
            mk_hit("r2", "A", 18, 25, '+'),
            mk_hit("r3", "B", 5, 15, '+'),
        ];
        let mut sorted = hits.clone();
        sorted.sort_by(|a, b| {
            (a.ref_name.as_str(), a.start, a.end).cmp(&(b.ref_name.as_str(), b.start, b.end))
        });
        let clusters = cluster_overlaps(sorted);
        assert_eq!(clusters.len(), 2); // A-cluster (2 hits) + B-cluster (1 hit)
        assert_eq!(clusters[0].len(), 2);
        assert_eq!(clusters[1].len(), 1);
    }

    #[test]
    fn test_build_read_maps_and_pairs() {
        // Two clusters on A and B, each with 3 hits (min support=3), sharing roots r1,r2,r3
        let mut hits: Vec<Hit> = Vec::new();
        for r in ["r1", "r2", "r3"] {
            hits.push(mk_hit(r, "A", 10, 20, '+'));
        }
        for r in ["r1", "r2", "r3"] {
            hits.push(mk_hit(r, "B", 30, 40, '+'));
        }
        hits.sort_by(|a, b| {
            (a.ref_name.as_str(), a.start, a.end).cmp(&(b.ref_name.as_str(), b.start, b.end))
        });
        let clusters = cluster_overlaps(hits);
        let read_maps = build_read_maps(clusters, 3);
        let p2r = build_pair_to_reads(&read_maps.read_to_regions);
        assert_eq!(p2r.len(), 1);
        let ((a, b), reads) = p2r.iter().next().unwrap();
        assert!(a.starts_with("A:") && b.starts_with("B:"));
        assert_eq!(reads.len(), 3);
        assert_eq!(read_maps.read_to_subreads.len(), 3);
    }
}
