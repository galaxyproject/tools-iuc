#!/usr/bin/env bash

fasta_file=$1
# Anolis_carolinensis.fa  miRNA   Anolis carolinensis
base_name=$(basename $fasta_file)
echo "$base_name miRNA Unknown species"
