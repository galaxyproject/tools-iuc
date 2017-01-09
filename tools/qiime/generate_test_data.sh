#!/usr/bin/env bash

# split_libraries_fastq
split_libraries_fastq.py \
    --sequence_read_fps 'test-data/forward_reads.fastq' \
    -o split_libraries \
    --mapping_fps 'test-data/map.tsv' \
    --barcode_read_fps 'test-data/barcodes.fastq' \
    --store_qual_scores \
    --store_demultiplexed_fastq \
    --max_bad_run_length 3 \
    --min_per_read_length_fraction 0.75 \
    --sequence_max_n 0 \
    --start_seq_id 0 \
    --barcode_type 'golay_12' \
    --max_barcode_errors 1.5
cp split_libraries/split_library_log.txt 'test-data/split_fastq_libraries_log.txt'
cp split_libraries/histograms.txt 'test-data/split_fastq_libraries_histograms.tabular'
cp split_libraries/seqs.fna 'test-data/split_fastq_libraries_sequences.fasta'
cp split_libraries/seqs.qual 'test-data/split_fastq_libraries_sequence_qualities.qual'
cp split_libraries/seqs.fastq 'test-data/split_fastq_libraries_demultiplexed_sequences.fastq'
rm -rf split_libraries
