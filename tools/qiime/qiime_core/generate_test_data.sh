#!/usr/bin/env bash

# validate_mapping_file
validate_mapping_file.py \
    -m 'test-data/validate_mapping_file/map.tsv' \
    -o validate_mapping_file_output \
    -c '_'
cp validate_mapping_file_output/*.html 'test-data/validate_mapping_file/map.tsv.html'
cp validate_mapping_file_output/*.log 'test-data/validate_mapping_file/map.tsv.log'
cp validate_mapping_file_output/*corrected.txt 'test-data/validate_mapping_file/map.tsv_corrected.txt'
rm -rf validate_mapping_file_output

# split_libraries_fastq
split_libraries_fastq.py \
    --sequence_read_fps 'test-data/split_libraries_fastq/forward_reads.fastq' \
    -o split_libraries \
    --mapping_fps 'test-data/map.tsv' \
    --barcode_read_fps 'test-data/split_libraries_fastq/barcodes.fastq' \
    --store_qual_scores \
    --store_demultiplexed_fastq \
    --max_bad_run_length 3 \
    --min_per_read_length_fraction 0.75 \
    --sequence_max_n 0 \
    --start_seq_id 0 \
    --barcode_type 'golay_12' \
    --max_barcode_errors 1.5
cp split_libraries/histograms.txt 'test-data/split_libraries_fastq/histograms.tabular'
cp split_libraries/seqs.fna 'test-data/split_libraries_fastq/sequences.fasta'
cp split_libraries/seqs.qual 'test-data/split_libraries_fastq/sequence_qualities.qual'
cp split_libraries/seqs.fastq 'test-data/split_libraries_fastq/demultiplexed_sequences.fastq'
rm -rf split_libraries

# split_libraries
split_libraries.py \
    --map 'test-data/split_libraries/mapping_file.txt' \
    -o split_libraries \
    --fasta 'test-data/split_libraries/reads_1.fna,test-data/split_libraries/reads_2.fna' \
    --qual 'test-data/split_libraries/reads_1.qual,test-data/split_libraries/reads_2.qual' \
    --min_qual_score 25 \
    --qual_score_window 0 \
    --record_qual_scores \
    --min_seq_length 200 \
    --max_seq_length 1000 \
    --max_ambig 6 \
    --max_homopolymer 6 \
    --max_primer_mismatch 0 \
    --barcode_type 'golay_12' \
    --max_barcode_errors 1.5 \
    --start_numbering_at 1
cp split_libraries/seqs.fna 'test-data/split_libraries/seqs.fna'
cp split_libraries/split_library_log.txt 'test-data/split_libraries/split_library_log'
cp split_libraries/histograms.txt 'test-data/split_libraries/histograms.txt'
cp split_libraries/seqs_filtered.qual 'test-data/split_libraries/seqs_filtered.qual'
rm -rf split_libraries

# pick_open_reference_otus
pick_open_reference_otus.py \
    --input_fps 'test-data/pick_open_reference_otus/sequences.fasta' \
    -o pick_open_reference_otus_1 \
    --reference_fp 'test-data/gg_13_8_79_otus.fasta' \
    --otu_picking_method 'uclust' \
    --new_ref_set_id 'New' \
    --parallel \
    --percent_subsample '0.001' \
    --prefilter_percent_id '0.0' \
    --minimum_failure_threshold '100000' \
    --min_otu_size '2'
cp pick_open_reference_otus_1/final_otu_map.txt 'test-data/pick_open_reference_otus/1_final_otu_map.txt'
cp pick_open_reference_otus_1/final_otu_map_mc*.txt 'test-data/pick_open_reference_otus/1_final_otu_map_mc.txt'
cp pick_open_reference_otus_1/rep_set.tre 'test-data/pick_open_reference_otus/1_rep_set_tree.tre'
rm -rf pick_open_reference_otus_1

pick_open_reference_otus.py \
    --input_fps 'test-data/pick_open_reference_otus/sequences.fasta' \
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
cp pick_open_reference_otus_2/final_otu_map.txt 'test-data/pick_open_reference_otus/2_final_otu_map.txt'
cp pick_open_reference_otus_2/final_otu_map_mc*.txt 'test-data/pick_open_reference_otus/2_final_otu_map_mc.txt'
rm -rf pick_open_reference_otus_2

pick_open_reference_otus.py \
    --input_fps 'test-data/pick_open_reference_otus/sequences.fasta' \
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
cp pick_open_reference_otus_3/final_otu_map.txt 'test-data/pick_open_reference_otus/3_final_otu_map.txt'
cp pick_open_reference_otus_3/final_otu_map_mc*.txt 'test-data/pick_open_reference_otus/3_final_otu_map_mc.txt'
cp pick_open_reference_otus_3/rep_set.tre 'test-data/pick_open_reference_otus/3_rep_set_tree.tre'
rm -rf pick_open_reference_otus_3

# core_diversity_analyses
# Data are from test data in https://github.com/biocore/qiime
core_diversity_analyses.py \
    --input_biom_fp 'test-data/core_diversity_analyses/otu_table.biom' \
    -o core_diversity_analyses_1 \
    --mapping_fp 'test-data/core_diversity_analyses/map.txt' \
    --sampling_depth 22 \
    --tree_fp 'test-data/core_diversity_analyses/rep_set.tre'
cp core_diversity_analyses_1/bdiv_even22/unweighted_unifrac_pc.txt 'test-data/core_diversity_analyses/unweighted_unifrac_pc.txt'
rm -rf core_diversity_analyses_1

core_diversity_analyses.py \
    --input_biom_fp 'test-data/core_diversity_analyses/otu_table.biom' \
    -o core_diversity_analyses_2 \
    --mapping_fp 'test-data/core_diversity_analyses/map.txt' \
    --sampling_depth 22 \
    --nonphylogenetic_diversity \
    --suppress_taxa_summary \
    --suppress_beta_diversity \
    --suppress_alpha_diversity \
    --suppress_group_significance
rm -rf core_diversity_analyses_2

# summarize_taxa
cp 'test-data/core_diversity_analyses/otu_table.biom' 'test-data/summarize_taxa/otu_table.biom'
cp 'test-data/core_diversity_analyses/map.txt' 'test-data/summarize_taxa/map.txt'

summarize_taxa.py \
    -i 'test-data/summarize_taxa/otu_table.biom' \
    -o summarize_taxa_1 \
    -L '2,3,4,5,6' \
    -m 'test-data/summarize_taxa/map.txt' \
    --md_identifier "taxonomy" \
    --delimiter ";"
cp summarize_taxa_1/*_L2.txt "test-data/summarize_taxa/1_L2.txt"
cp summarize_taxa_1/*_L3.txt "test-data/summarize_taxa/1_L3.txt"
cp summarize_taxa_1/*_L4.txt "test-data/summarize_taxa/1_L4.txt"
cp summarize_taxa_1/*_L5.txt "test-data/summarize_taxa/1_L5.txt"
cp summarize_taxa_1/*_L6.txt "test-data/summarize_taxa/1_L6.txt"
rm -rf summarize_taxa_1

summarize_taxa.py \
    -i 'test-data/summarize_taxa/otu_table.biom' \
    -o summarize_taxa_2 \
    -L '3,6' \
    --md_identifier "taxonomy" \
    --delimiter ";"
cp summarize_taxa_2/*_L3.txt "test-data/summarize_taxa/2_L3.txt"
cp summarize_taxa_2/*_L6.txt "test-data/summarize_taxa/2_L6.txt"
rm -rf summarize_taxa_2

# make_emperor
cp 'test-data/core_diversity_analyses/unweighted_unifrac_pc.txt' 'test-data/make_emperor/unweighted_unifrac_pc.txt'
cp 'test-data/core_diversity_analyses/map.txt' 'test-data/make_emperor/map.txt'
cp 'test-data/summarize_taxa/2_L3.txt' 'test-data/make_emperor/2_L3.txt'

make_emperor.py \
    --input_coords 'test-data/make_emperor/unweighted_unifrac_pc.txt' \
    -o make_emperor_1 \
    --map_fp 'test-data/make_emperor/map.txt' \
    --number_of_axes '10' \
    --add_unique_columns \
    --number_of_segments 8
rm -rf make_emperor_1

make_emperor.py \
    --input_coords 'test-data/make_emperor/unweighted_unifrac_pc.txt' \
    -o make_emperor_2 \
    --map_fp 'test-data/make_emperor/map.txt' \
    --number_of_axes '10' \
    --add_unique_columns \
    --number_of_segments 8 \
    --taxa_fp 'test-data/make_emperor/2_L3.txt' \
    --n_taxa_to_keep 10
rm -rf make_emperor_2

#alpha_rarefaction
alpha_rarefaction.py \
    --otu_table_fp "test-data/alpha_rarefaction/otu_table.biom" \
    --mapping_fp "test-data/alpha_rarefaction/mapping_file.txt" \
    -o alpha_rarefaction \
    --num_steps '2' \
    --tree_fp "test-data/alpha_rarefaction/rep_set.tre" \
    --min_rare_depth '10' \
    --max_rare_depth '50' \
    --retain_intermediate_files
rm -rf alpha_rarefaction

##beta_diversity
beta_diversity.py \
    --input_path 'test-data/beta_diversity/otu_table.biom' \
    -o beta_diversity_1 \
    --metrics 'unweighted_unifrac,weighted_unifrac' \
    --tree_path 'test-data/beta_diversity/rep_set.tre'
md5 'beta_diversity_1/unweighted_unifrac_otu_table.txt'
md5 'beta_diversity_1/weighted_unifrac_otu_table.txt'
rm -rf beta_diversity_1

beta_diversity.py \
    --input_path 'test-data/beta_diversity/otu_table.biom' \
    -o beta_diversity_2 \
    --metrics 'abund_jaccard,binary_chisq,binary_chord,binary_euclidean,binary_hamming,binary_jaccard,binary_lennon,binary_ochiai,binary_pearson,binary_sorensen_dice,bray_curtis,canberra,chisq,chord,euclidean,gower,hellinger,kulczynski,manhattan,morisita_horn,pearson,soergel,spearman_approx,specprof,unifrac_g,unifrac_g_full_tree,unweighted_unifrac,unweighted_unifrac_full_tree,weighted_normalized_unifrac,weighted_unifrac' \
    --tree_path 'test-data/beta_diversity/rep_set.tre'
md5 'beta_diversity_2/canberra_otu_table.txt'
md5 'beta_diversity_2/pearson_otu_table.txt'
rm -rf beta_diversity_2

# compare_categories
compare_categories.py \
    --method 'adonis' \
    --input_dm 'test-data/compare_categories/unweighted_unifrac_dm.txt' \
    --mapping_file 'test-data/compare_categories/map.txt' \
    --categories 'Treatment' \
    -o compare_categories_1 \
    --num_permutations '999'
cp compare_categories_1/adonis_results.txt "test-data/compare_categories/adonis_results.txt"
rm -rf compare_categories_1

compare_categories.py \
    --method 'dbrda' \
    --input_dm 'test-data/compare_categories/unweighted_unifrac_dm.txt' \
    --mapping_file 'test-data/compare_categories/map.txt' \
    --categories 'Treatment' \
    -o compare_categories_2 \
    --num_permutations '99'
cp compare_categories_2/* "test-data/compare_categories/"
rm -rf compare_categories_2
