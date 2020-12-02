#!/bin/bash

# This script produces a small kraken2 database containing only a ~1kb portion each of a salmonella and ecoli genome
# It requires kraken2, art and entrez-direct (all available on bioconda)
kraken2-build --db test_db --download_taxonomy
mv test_db/taxonomy/nucl_gb.accession2taxid test_db/taxonomy/nucl_gb.accession2taxid_full
grep -e 'NC_003198.1' -e 'NC_011750.1' test_db/taxonomy/nucl_gb.accession2taxid_full > test_db/taxonomy/nucl_gb.accession2taxid
esearch -db nucleotide -query "NC_003198.1" | efetch -format fasta > NC_003198.1.fasta
esearch -db nucleotide -query "NC_011750.1" | efetch -format fasta > NC_011750.1.fasta
head -n 14 NC_003198.1.fasta > NC_003198.1_1kb.fasta
head -n 14 NC_011750.1.fasta > NC_011750.1_1kb.fasta
kraken2-build --db test_db --add-to-library NC_003198.1_1kb.fasta
kraken2-build --db test_db --add-to-library NC_011750.1_1kb.fasta
kraken2-build --db test_db --build

# Simulate 100bp reads from ~1kb portions of genomes
art_illumina -sam -i NC_011750.1_1kb.fasta -p -m 300 -f 10 -s 10 -l 100 -o NC_011750.1_simulated_R
art_illumina -sam -i NC_003198.1_1kb.fasta -p -m 300 -f 10 -s 10 -l 100 -o NC_003198.1_simulated_R

# Generate kraken reports
kraken2 --db test_db --report NC_011750.1_simulated_kraken_report.txt --paired NC_011750.1_simulated_R1.fastq NC_011750.1_simulated_R2.fastq
kraken2 --db test_db --report NC_003198.1_simulated_kraken_report.txt --paired NC_003198.1_simulated_R1.fastq NC_003198.1_simulated_R2.fastq

# Build bracken kmer distribution files using default kmer-len=35 and read-len=100
bracken-build -d test_db

# 
# est_abundance.py --kmer_distr test_db/database100mers.kmer_distrib --level S -i NC_003198.1_simulated_kraken_report.txt -o NC_003198.1_simulated_bracken_report.txt
# est_abundance.py --kmer_distr test_db/database100mers.kmer_distrib --level S -i NC_011750.1_simulated_kraken_report.txt -o NC_011750.1_simulated_bracken_report.txt
