#!/bin/bash

# 1. Create a database of unique kmers from the parental genomes
cd examples
mkdir -p ref f1_reads f1_classify parent_reads parent_classify
faSplit byname parents.genome.fa ref/
klassify build ref/*.fa -o kmers.bc

# 2. Classify the progeny (e.g. F1) reads based on the unique kmers
faSplit about f1_reads.fa 2000000000 f1_reads/
klassify classify kmers.bc f1_reads/*.fa -o f1_classify

# 3. Map ‘chimeric’ progeny reads to the parents reference
klassify extract f1_classify.filtered.tsv f1_reads/*.fa -o f1_classify.fa
minimap2 -t 80 -ax map-hifi --eqx --secondary=no parents.genome.fa f1_classify.fa \
    --split-prefix f1_classify | samtools sort -@ 8 -o f1_classify.bam

# 4. Repeat the steps using the parental reads
faSplit about parent_reads.fa 2000000000 parent_reads/
klassify classify kmers.bc parent_reads/*.fa -o parent_classify
klassify extract parent_classify.filtered.tsv parent_reads/*.fa -o parent_classify.fa
minimap2 -t 80 -ax map-hifi --eqx --secondary=no parents.genome.fa parent_classify.fa \
    --split-prefix parent_classify | samtools sort -@ 8 -o parent_classify.bam

# 5. Using parent reads as ‘control’, identify the ‘chimeric’ regions that show
#    up with F1 reads, but NOT with parent reads (so we are not affected by
#    assembly errors)
klassify regions f1_classify.bam parent_classify.bam
