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

1. Create a database of unique kmers from the parental genomes

```console
klassify build ref/*.fa
```

This generates a file called `singleton_kmers.bc` that indexes all the unique kmers.

2. Classify the progeny (e.g. F1) reads based on the unique kmers

```console
klassify classify singleton_kmers.bc f1_reads/*.fa
```

3. Map ‘chimeric’ progeny reads to the reference

```console
minimap2 -t 80 -ax map-hifi --eqx --secondary=no ref/parents.genome.fasta \
    | samtools sort -@ 8 -o chimeric_f1.parents.bam
```

4. Repeat the steps using the parental reads

```console
klassify classify singleton_kmers.bc parent_reads/*.fa
minimap2 -t 80 -ax map-hifi --eqx --secondary=no ref/parents.genome.fasta \
    | samtools sort -@ 8 -o chimeric_parents.parents.bam
```

5. Using parent reads as ‘control’, identify the ‘chimeric’ regions that show up with F1 reads, but NOT with parent reads (so we are not affected by assembly errors)

```console
klassify compare chimeric_f1.parents.bam chimeric_parents.parents.bam
```
