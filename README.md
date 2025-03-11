# klassify

[![Crates.io](https://img.shields.io/crates/v/klassify.svg)](https://crates.io/crates/klassify)
[![Github Actions](https://github.com/tanghaibao/klassify/actions/workflows/rust.yml/badge.svg)](https://github.com/tanghaibao/klassify/actions)

![klassify-logo](https://www.dropbox.com/scl/fi/bjvfamep0aoxka0dcg2zi/klassify-logo.png?rlkey=8vmvacehs2amuaoi0gvgyh28r&st=ohygf458&raw=1)

Classify chimeric reads based on unique k-mer contents and identify the
breakpoint locations.

The breakpoints can be due to:

- Recombination / crossover events
- Structural variations

While there are many tools that can identify structural variations, this tool
is designed to compare progeny (e.g. F1) reads to the parental genome. The key
idea is an extension to the trio binning approach, where we use the unique kmers
from each chromosome/contig of the parental genomes to classify the reads that
bridge two different chromosomes/contigs.

Following are examples of recominant reads identified by this tool:

![recombinant-read](https://www.dropbox.com/scl/fi/tduxwsh0wcy2zdw8zopdm/recombinant-reads.png?rlkey=xci43gwwy84dbcdvs2n7ekk18&st=sewwc9s0&raw=1)

## Installation

If you don't have [Rust](https://rustup.rs/) installed:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

With Rust installed, you can just install the software with `cargo`.

```bash
cargo install klassify
```

Typical install time on a desktop computer is ~1 minute.

Additional dependencies include:

- [`minimap2`](https://github.com/lh3/minimap2)
- [`samtools`](https://github.com/samtools/samtools)
- [`faSplit`](https://hgdownload.soe.ucsc.edu/admin/exe/)
- [`mosdepth`](https://github.com/brentp/mosdepth)

## Supported Operating systems

We have tested latest version (`0.1.4` and above) on the following OS:

- Linux (Redhat Enterprise 7.6)

## Usage

Suppose you have 3 input files, with a toy example available in `examples`:

- `parents.genome.fa`: the parental genomes
- `f1_reads.fa`: the progeny reads
- `parent_reads.fa`: the parental reads

1. Create a database of unique kmers from the parental genomes

```console
cd examples
mkdir -p ref f1_reads f1_classify parent_reads parent_classify
faSplit byname parents.genome.fa ref/
klassify build ref/*.fa -o kmers.bc
```

This generates an index for all the unique kmers (present in a single contig/chromosome).

2. Classify the progeny (e.g. F1) reads based on the unique kmers

```console
faSplit about f1_reads.fa 2000000000 f1_reads/
klassify classify kmers.bc f1_reads/*.fa -o f1_classify
```

3. Map ‘chimeric’ progeny reads to the parents reference

```console
klassify extract f1_classify.filtered.tsv f1_reads/*.fa -o f1_classify.fa
minimap2 -t 80 -ax map-hifi --eqx --secondary=no parents.genome.fa f1_classify.fa \
    --split-prefix f1_classify | samtools sort -@ 8 -o f1_classify.bam
```

4. Repeat the steps using the parental reads

```console
faSplit about parent_reads.fa 2000000000 parent_reads/
klassify classify kmers.bc parent_reads/*.fa -o parent_classify
klassify extract parent_classify.filtered.tsv parent_reads/*.fa -o parent_classify.fa
minimap2 -t 80 -ax map-hifi --eqx --secondary=no parents.genome.fa parent_classify.fa \
    --split-prefix parent_classify | samtools sort -@ 8 -o parent_classify.bam
```

5. Using parent reads as ‘control’, identify the ‘chimeric’ regions that show up with F1 reads, but NOT with parent reads (so we are not affected by assembly errors)

```console
klassify regions f1_classify.bam parent_classify.bam
```

That's it! The breakpoint locations in the parental genomes are in
`f1_classify.mosdepth.regions.bed.regions.tsv`, where column 2 shows the supported
depth within each consecutive 10kb bin around the breakpoint (by default: at
least 5 supported reads):

```console
SoChr01B:70000-90000	12,5
SoChr01F:80000-90000	11
```

The breakpoint locations can then be visualized in IGV for read evidence in
`f1_classify.bam`, using `parents.genome.fa` as the reference.

Total expected run time on a desktop computer is ~1 minute.

## Algorithm

The KLASSIFY pipeline identifies the breakpoints using the set of F1 reads,
with parent reads as control. The breakpoints identified from the F1 reads were
then mapped back to the parent reference sequences to obtain precise coordinates.

The KLASSIFY algorithm works as follows:

1. Find unique k-mers that belong to each chromosome, e.g. SoChr01A, SoChr01B,
   etc.

2. Identify ‘chimeric’ F1 reads that contain unique k-mers that belong to at
   least 2 chromosomes (default: ≧300 unique k-mers on the read, A unique + B unique
   ≧50% of unique k-mers on the read, and B unique ≧10%)

3. Repeat step 2 similarly with parent reads

4. Using parent reads as ‘control’, identify the ‘chimeric’ regions that show up
   with at least 5 F1 reads, but not with parent reads (therefore unaffected by
   assembly errors or repeats)

5. Collect all ‘chimeric’ reads identified so far and split them into 2 parts.
   The reads are split by identifying the switch from one chromosome to another
   based on unique k-mers

6. Map the split reads to the reference sequences to identify parent regions
   where each part of the ‘chimeric’ reads separately map to

7. Pair the separate regions up to compile a candidate list of paired breakpoints

8. Use Integrated Genome Viewer (IGV) to proof the paired breakpoints. Label the
   breakpoint as either “Type I”, “Type II”, or “bad” (see next section for
   definition of types)