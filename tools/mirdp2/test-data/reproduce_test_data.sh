#!/bin/bash

# This script produces a small bowtie index
# It requires bowtie
gunzip genome.fasta.gz
mkdir test_db/tair10_index
bowtie-build -f genome.fasta test_db/tair10_index
