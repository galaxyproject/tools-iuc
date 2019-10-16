#!/bin/bash

# This script produces a small kraken2 database containing only a ~1kb portion each of a salmonella and ecoli genome
# It requires kraken2, and entrez-direct (available on bioconda)
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
