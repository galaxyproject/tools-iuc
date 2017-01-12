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

# pick_open_reference_otus
pick_open_reference_otus.py \
    --input_fps 'test-data/split_fastq_libraries_sequences.fasta' \
    -o pick_open_reference_otus_1 \
    --reference_fp 'test-data/gg_13_8_79_otus.fasta' \
    --otu_picking_method 'uclust' \
    --new_ref_set_id 'New' \
    --parallel \
    --percent_subsample '0.001' \
    --prefilter_percent_id '0.0' \
    --minimum_failure_threshold '100000' \
    --min_otu_size '2'
cp pick_open_reference_otus_1/final_otu_map.txt 'test-data/pick_open_reference_otus_1_final_otu_map.txt'
cp pick_open_reference_otus_1/final_otu_map_mc*.txt 'test-data/pick_open_reference_otus_1_final_otu_map_mc.txt'
cp pick_open_reference_otus_1/rep_set.fna 'test-data/pick_open_reference_otus_1_rep_set.fna'
cp pick_open_reference_otus_1/rep_set.tre 'test-data/pick_open_reference_otus_1_rep_set_tree.tre'
rm -rf pick_open_reference_otus_1

pick_open_reference_otus.py \
    --input_fps 'test-data/split_fastq_libraries_sequences.fasta' \
    -o pick_open_reference_otus_2 \
    --reference_fp 'test-data/gg_13_8_79_otus.fasta' \
    --otu_picking_method 'uclust' \
    --new_ref_set_id 'New' \
    --parallel \
    --percent_subsample '0.001' \
    --prefilter_percent_id '0.0' \
    --minimum_failure_threshold '100000' \
    --min_otu_size '3' \
    --suppress_taxonomy_assignment \
    --suppress_align_and_tree
cp pick_open_reference_otus_2/final_otu_map.txt 'test-data/pick_open_reference_otus_2_final_otu_map.txt'
cp pick_open_reference_otus_2/final_otu_map_mc*.txt 'test-data/pick_open_reference_otus_2_final_otu_map_mc.txt'
rm -rf pick_open_reference_otus_2

pick_open_reference_otus.py \
    --input_fps 'test-data/split_fastq_libraries_sequences.fasta' \
    -o pick_open_reference_otus_3 \
    --reference_fp 'test-data/gg_13_8_79_otus.fasta' \
    --otu_picking_method 'uclust' \
    --new_ref_set_id 'New' \
    --parallel \
    --percent_subsample '0.001' \
    --prefilter_percent_id '0.0' \
    --minimum_failure_threshold '100000' \
    --min_otu_size '10' \
    --suppress_taxonomy_assignment
cp pick_open_reference_otus_3/final_otu_map.txt 'test-data/pick_open_reference_otus_3_final_otu_map.txt'
cp pick_open_reference_otus_3/final_otu_map_mc*.txt 'test-data/pick_open_reference_otus_3_final_otu_map_mc.txt'
cp pick_open_reference_otus_3/rep_set.tre 'test-data/pick_open_reference_otus_3_rep_set_tree.tre'
rm -rf pick_open_reference_otus_3

# core_diversity_analyses
# Data are from test data in https://github.com/biocore/qiime
core_diversity_analyses.py \
    --input_biom_fp 'test-data/core_diversity_analyses_otu_table.biom' \
    -o core_diversity_analyses_1 \
    --mapping_fp 'test-data/core_diversity_analyses_map.txt' \
    --sampling_depth 22 \
    --tree_fp 'test-data/core_diversity_analyses_rep_set.tre'
rm -rf core_diversity_analyses_1

core_diversity_analyses.py \
    --input_biom_fp 'test-data/core_diversity_analyses_otu_table.biom' \
    -o core_diversity_analyses_2 \
    --mapping_fp 'test-data/core_diversity_analyses_map.txt' \
    --sampling_depth 22 \
    --nonphylogenetic_diversity \
    --suppress_taxa_summary \
    --suppress_beta_diversity \
    --suppress_alpha_diversity \
    --suppress_group_significance
rm -rf core_diversity_analyses_2
