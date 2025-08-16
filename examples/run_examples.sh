#!/bin/bash

# 1. Create a database of unique kmers from the parental genomes
cd examples || exit 1
mkdir -p f1_classify parent_classify
klassify build parents.genome.fa -o kmers.bc

# 2. Classify the progeny (e.g. F1) reads based on the unique kmers, extract the
#    reads that are classified as ‘chimeric’ and map them to the parents reference
klassify classify kmers.bc f1_reads.fa -o f1_classify
klassify extract f1_classify.filtered.tsv f1_reads.fa -o f1_classify.fa
minimap2 -t 80 -ax map-hifi --eqx --secondary=no parents.genome.fa f1_classify.fa \
    | samtools sort -@ 8 -o f1_classify.bam

# 3. Repeat the steps using the parental reads
klassify classify kmers.bc parent_reads.fa -o parent_classify
klassify extract parent_classify.filtered.tsv parent_reads.fa -o parent_classify.fa
minimap2 -t 80 -ax map-hifi --eqx --secondary=no parents.genome.fa parent_classify.fa \
    | samtools sort -@ 8 -o parent_classify.bam

# 4. Using parent reads as ‘control’, identify the ‘chimeric’ regions that show
#    up with F1 reads, but NOT with parent reads (so we are not affected by
#    assembly errors)
klassify regions f1_classify.bam parent_classify.bam
