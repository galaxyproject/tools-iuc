#!/usr/bin/bash

# E. coli locus b0842 (b0842.fasta.gz) downloaded from Enterobase E. coli cgMLST scheme
# requires: wget, kma, bwa, samtools, bedtools

gunzip b0842.fasta.gz

# Take first 5 alleles to reduce size of test data
mkdir ecoli_cgMLST
head -n 10 b0842.fasta > ecoli_cgMLST/ecoli_b0842_1to5.fasta

kma index -k 8 -i ecoli_cgMLST/ecoli_b0842_1to5.fasta -o ecoli_cgMLST/ecoli_b0842_1to5

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR884/ERR884056/ERR884056_1.fastq.gz

# Use bwa to map reads to reduced E. coli locus b0842
# and extract only mapped reads (to reduce size of test dataset)
bwa index ecoli_cgMLST/ecoli_b0842_1to5.fasta

bwa mem ecoli_cgMLST/ecoli_b0842_1to5.fasta ERR884056_1.fastq.gz -o ERR884056_1_ecoli_b0842_1to5.sam

samtools view ERR884056_1_ecoli_b0842_1to5.sam -bo ERR884056_1_ecoli_b0842_1to5.bam

# Select mapped reads
samtools view -b -F 4 ERR884056_1_ecoli_b0842_1to5.bam > ERR884056_1_ecoli_b0842_1to5.mapped.bam

samtools sort -n ERR884056_1_ecoli_b0842_1to5.mapped.bam -o ERR884056_1_ecoli_b0842_1to5.mapped.sort.bam

bedtools bamtofastq -i ERR884056_1_ecoli_b0842_1to5.mapped.sort.bam -fq ERR884056_ecoli_b0842.mapped_R1.fastq

