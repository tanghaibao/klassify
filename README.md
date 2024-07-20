# klassify

![klassify-logo](https://www.dropbox.com/scl/fi/bjvfamep0aoxka0dcg2zi/klassify-logo.png?rlkey=8vmvacehs2amuaoi0gvgyh28r&st=ohygf458&raw=1)

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

Following are examples of recominant reads identified by this tool:

![recombinant-read](https://www.dropbox.com/scl/fi/tduxwsh0wcy2zdw8zopdm/recombinant-reads.png?rlkey=xci43gwwy84dbcdvs2n7ekk18&st=sewwc9s0&raw=1)

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
mkdir ref
faSplit byname parents.genome.fa ref/
klassify build ref/*.fa
```

This generates a file called `singleton_kmers.bc` that indexes all the unique kmers.

2. Classify the progeny (e.g. F1) reads based on the unique kmers

```console
mkdir f1_reads f1_classify
faSplit about f1_reads.fa 2000000000 f1_reads/
klassify classify singleton_kmers.bc f1_reads/*.fa
mv *.read_classfications.tsv f1_classify/
python scripts/classify_reads.py f1_classify
```

3. Map ‘chimeric’ progeny reads to the parents reference

```console
klassify extract f1_classify.filtered.tsv f1_reads/*.fa
cat *.extracted.fa > f1_classify.fa
rm *.extracted.fa
minimap2 -t 80 -ax map-hifi --eqx --secondary=no ref/parents.genome.fa f1_classify.fa \
    --split-prefix f1_classify | samtools sort -@ 8 -o f1_classify.bam
```

4. Repeat the steps using the parental reads

```console
mkdir parent_reads parent_classify
faSplit about parent_reads.fa 1000000000 parent_reads/
klassify classify singleton_kmers.bc parent_reads/*.fa
mv *.read_classfications.tsv parent_classify/
python scripts/classify_reads.py parent_classify
klassify extract parent_classify.filtered.tsv parent_reads/*.fa
cat *.extracted.fa > parent_classify.fa
rm *.extracted.fa
minimap2 -t 80 -ax map-hifi --eqx --secondary=no ref/parents.genome.fa parent_classify.fa \
    --split-prefix parent_classify | samtools sort -@ 8 -o parent_classify.bam
```

5. Using parent reads as ‘control’, identify the ‘chimeric’ regions that show up with F1 reads, but NOT with parent reads (so we are not affected by assembly errors)

```console
samtools index f1_classify.bam
samtools index parent_classify.bam
mosdepth -t 8 -n --by 10000 f1_classify.mosdepth f1_classify.bam
mosdepth -t 8 -n --by 10000 parent_classify.mosdeth parent_classify.bam
python scripts/generate_regions.py f1_classify.bam parent_classify.bam
```
