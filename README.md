# klassify

Classify chimeric reads based on unique kmer contents and identify the
breakpoint locations.

The breakpoints can be due to:

- Recombination / crossover events
- Structrual variations

While there are many tools that can identify structural variations, this tool
is designed to compare progeny (e.g. F1) reads to the parental genome. The key
idea is an extension to the trio binning approach, where we use the unique kmers
from each chromosome/contig of the parental genomes to classify the reads that
bridge two different chromosomes/contigs.

## Installation

```bash
cargo build -r
```

## Usage

Suppose you have 3 input files:

- `parents.genome.fa`: the parental genomes
- `f1_reads.fa`: the progeny reads
- `parent_reads.fa`: the parental reads

1. Create a database of unique kmers from the parental genomes

```console
mkdir ref/
faSplit byname parents.genome.fa ref/
klassify build ref/*.fa
```

This generates a file called `singleton_kmers.bc` that indexes all the unique kmers.

2. Classify the progeny (e.g. F1) reads based on the unique kmers

```console
mkdir f1_reads
faSplit about f1_reads.fa 1000000000 f1_reads/
klassify classify singleton_kmers.bc f1_reads/*.fa
python classify_reads.py
```

3. Map ‘chimeric’ progeny reads to the reference

```console
klassify extract chimeric_f1_reads.tsv f1_reads/*.fa
cat *.extracted.fa > chimeric_f1_reads.fa
minimap2 -t 80 -ax map-hifi --eqx --secondary=no ref/parents.genome.fa chimeric_f1_reads.fa \
    | samtools sort -@ 8 -o chimeric_f1_reads.parents.bam
```

4. Repeat the steps using the parental reads

```console
mkdir parent_reads
faSplit about parent_reads.fa 1000000000 parent_reads/
klassify classify singleton_kmers.bc parent_reads/*.fa
python classify_reads.py
klassify extract chimeric_parent_reads.tsv parent_reads/*.fa
cat *.extracted.fa > chimeric_parent_reads.fa
minimap2 -t 80 -ax map-hifi --eqx --secondary=no ref/parents.genome.fa chimeric_parent_reads.fa \
    | samtools sort -@ 8 -o chimeric_parent_reads.parents.bam
```

5. Using parent reads as ‘control’, identify the ‘chimeric’ regions that show up with F1 reads, but NOT with parent reads (so we are not affected by assembly errors)

```console
samtools index chimeric_f1_reads.parents.bam
samtools index chimeric_parent_reads.parents.bam
mosdepth -t 8 -n --by 10000 chimeric_f1_reads.mosdeth chimeric_f1_reads.parents.bam
mosdepth -t 8 -n --by 10000 chimeric_parent_reads.mosdeth chimeric_parent_reads.parents.bam
```
