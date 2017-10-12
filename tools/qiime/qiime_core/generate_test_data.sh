#!/usr/bin/env bash

# Data are from test data in https://github.com/biocore/qiime

# align_seqs
align_seqs.py \
    --input_fasta_fp 'test-data/align_seqs/unaligned.fna' \
    -o 'align_seqs_pynast_uclust' \
    --alignment_method 'pynast' \
    --pairwise_alignment_method 'uclust' \
    --template_fp 'test-data/align_seqs/core_set_aligned.fasta.imputed' \
    --min_percent_id '0.75'

align_seqs.py \
    --input_fasta_fp 'test-data/align_seqs/unaligned.fna' \
    -o 'align_seqs_pynast_muscle' \
    --alignment_method 'pynast' \
    --pairwise_alignment_method 'muscle' \
    --min_length '50' \
    --min_percent_id '0.75'

align_seqs.py \
    --input_fasta_fp 'test-data/align_seqs/unaligned.fna' \
    -o 'align_seqs_pynast_pair_hmm' \
    --alignment_method 'pynast' \
    --pairwise_alignment_method 'pair_hmm' \
    --min_percent_id '0.75'

#align_seqs.py \
#    --input_fasta_fp 'test-data/align_seqs/unaligned.fna' \
#    -o 'align_seqs_pynast_clustal' \
#    --alignment_method 'pynast' \
#    --pairwise_alignment_method 'clustal' \
#    --min_percent_id '0.75'

align_seqs.py \
    --input_fasta_fp 'test-data/align_seqs/unaligned.fna' \
    -o 'align_seqs_pynast_blast' \
    --alignment_method 'pynast' \
    --pairwise_alignment_method 'blast' \
    --min_percent_id '0.75'

align_seqs.py \
    --input_fasta_fp 'test-data/align_seqs/unaligned.fna' \
    -o 'align_seqs_pynast_mafft' \
    --alignment_method 'pynast' \
    --pairwise_alignment_method 'mafft' \
    --min_percent_id '0.75'

#align_seqs.py \
#    --input_fasta_fp 'test-data/align_seqs/unaligned.fna' \
#    -o 'align_seqs_infernal' \
#    --alignment_method 'infernal' \
#    --template_fp 'test-data/align_seqs/seed.16s.reference_model.sto' \
#    --min_percent_id '0.75'

#align_seqs.py \
#    --input_fasta_fp 'test-data/align_seqs/unaligned.fna' \
#    -o 'align_seqs_clustalw' \
#    --alignment_method 'clustalw' \
#    --min_percent_id '0.75'

align_seqs.py \
    --input_fasta_fp 'test-data/align_seqs/unaligned.fna' \
    -o 'align_seqs_muscle' \
    --alignment_method 'muscle' \
    --min_percent_id '0.75'

align_seqs.py \
    --input_fasta_fp 'test-data/align_seqs/unaligned.fna' \
    -o 'align_seqs_mafft' \
    --alignment_method 'mafft' \
    --min_percent_id '0.75'

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

# assign_taxonomy
assign_taxonomy.py \
    --input_fasta_fp 'test-data/assign_taxonomy/uclust_input_seqs.fasta' \
    --assignment_method 'uclust' \
    --min_consensus_fraction '0.51' \
    --similarity '0.9' \
    --uclust_max_accepts '3' \
    -o assign_taxonomy_uclust
ls assign_taxonomy_uclust
md5sum 'assign_taxonomy_uclust/uclust_input_seqs_tax_assignments.txt'
rm -rf assign_taxonomy_uclust

assign_taxonomy.py \
    --input_fasta_fp 'test-data/assign_taxonomy/mothur_repr_set_seqs.fasta' \
    --id_to_taxonomy_fp 'test-data/assign_taxonomy/mothur_id_to_taxonomy.txt' \
    --assignment_method 'mothur' \
    --reference_seqs_fp 'test-data/assign_taxonomy/mothur_ref_seq_set.fna' \
    --confidence '0.5' \
    -o assign_taxonomy_mothur
ls assign_taxonomy_mothur
md5sum 'assign_taxonomy_mothur/mothur_repr_set_seqs_tax_assignments.txt'
rm -rf assign_taxonomy_mothur

assign_taxonomy.py \
    --input_fasta_fp 'test-data/assign_taxonomy/mothur_repr_set_seqs.fasta' \
    --id_to_taxonomy_fp 'test-data/assign_taxonomy/mothur_id_to_taxonomy.txt' \
    --assignment_method 'mothur' \
    --reference_seqs_fp 'test-data/assign_taxonomy/mothur_ref_seq_set.fna' \
    --blast_e_value '0.001' \
    -o assign_taxonomy_blast
ls assign_taxonomy_blast
md5sum 'assign_taxonomy_blast/mothur_repr_set_seqs_tax_assignments.txt'
rm -rf assign_taxonomy_blast

#assign_taxonomy.py \
#    --input_fasta_fp 'test-data/assign_taxonomy/rdp_input_seqs.fasta' \
#    --id_to_taxonomy_fp 'test-data/assign_taxonomy/rdp_id_to_taxonomy.txt' \
#    --assignment_method 'rdp' \
#    --confidence '3' \
#    -o assign_taxonomy_rdp

#assign_taxonomy.py \
#    --input_fasta_fp 'test-data/assign_taxonomy/rtax_ref_seq_set.fna' \
#    --id_to_taxonomy_fp 'test-data/assign_taxonomy/rtax_id_to_taxonomy.txt' \
#    --assignment_method 'rtax' \
#    --read_1_seqs_fp 'test-data/assign_taxonomy/read_1.seqs.fna' \
#    --read_2_seqs_fp 'test-data/assign_taxonomy/read_2.seqs.fna'  \
#    --single_ok \
#    --no_single_ok_generic \
#    --read_id_regex "\S+\s+(\S+)" \
#    --amplicon_id_regex "(\S+)\s+(\S+?)\/" \
#    --header_id_regex "\S+\s+(\S+?)\/" \
#    -o assign_taxonomy_rtax
#ls assign_taxonomy_rtax

assign_taxonomy.py \
    --input_fasta_fp 'test-data/assign_taxonomy/mothur_ref_seq_set.fna' \
    --assignment_method 'sortmerna' \
    --min_consensus_fraction "0.51" \
    --similarity "0.9" \
    --sortmerna_e_value "1.0" \
    --sortmerna_coverage "0.9" \
    --sortmerna_best_N_alignments "5" \
    -o assign_taxonomy_sortmerna
ls assign_taxonomy_sortmerna
md5sum 'assign_taxonomy_sortmerna/mothur_ref_seq_set_tax_assignments.txt'
md5sum 'assign_taxonomy_sortmerna/sortmerna_map.blast'
rm -rf assign_taxonomy_sortmerna

#beta_diversity
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

#beta_diversity_through_plots
beta_diversity_through_plots.py \
    --otu_table_fp 'test-data/beta_diversity_through_plots/otu_table.biom' \
    --mapping_fp 'test-data/beta_diversity_through_plots/map.txt' \
    --output_dir beta_diversity_through_plots \
    --tree_fp 'test-data/beta_diversity_through_plots/rep_set.tre' \
    --parallel
cp beta_diversity_through_plots/unweighted_unifrac_dm.txt 'test-data/beta_diversity_through_plots/'
cp beta_diversity_through_plots/unweighted_unifrac_pc.txt 'test-data/beta_diversity_through_plots/'
cp beta_diversity_through_plots/weighted_unifrac_dm.txt 'test-data/beta_diversity_through_plots/'
cp beta_diversity_through_plots/weighted_unifrac_pc.txt 'test-data/beta_diversity_through_plots/'
rm -rf beta_diversity_through_plots

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

# core_diversity_analyses
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

# extract_barcodes
extract_barcodes.py \
    --fastq1 'test-data/extract_barcodes/inseqs.fastq' \
    --input_type 'barcode_single_end' \
    -o extract_barcodes_1 \
    --bc1_len '6' \
    --rev_comp_bc1
rm -rf extract_barcodes_1

extract_barcodes.py \
    --fastq1 'test-data/extract_barcodes/inseqs_R1.fastq' \
    --input_type 'barcode_paired_end' \
    --fastq2 'test-data/extract_barcodes/inseqs_R2.fastq' \
    -o extract_barcodes_2 \
    --bc1_len '6' \
    --bc2_len '6'
rm -rf extract_barcodes_2

extract_barcodes.py \
    --fastq1 'test-data/extract_barcodes/inseqs_R1.fastq' \
    --input_type 'barcode_paired_end' \
    --fastq2 'test-data/extract_barcodes/inseqs_R2.fastq' \
    -o extract_barcodes_3 \
    --bc1_len '6' \
    --bc2_len '6' \
    --mapping_fp 'test-data/extract_barcodes/mapping_data.txt' \
    --attempt_read_reorientation \
    --disable_header_match
rm -rf extract_barcodes_3

extract_barcodes.py \
    --fastq1 'test-data/extract_barcodes/inseqs_R1.fastq' \
    --input_type 'barcode_paired_stitched' \
    -o extract_barcodes_4 \
    --bc1_len '6' \
    --bc2_len '8' \
    --rev_comp_bc1 \
    --rev_comp_bc2
rm -rf extract_barcodes_4

extract_barcodes.py \
    --fastq1 'test-data/extract_barcodes/inseqs_R1.fastq' \
    --input_type 'barcode_in_label' \
    --char_delineator '#' \
    -o extract_barcodes_5 \
    --bc1_len '6'
rm -rf extract_barcodes_5

# filter_alignment
filter_alignment.py \
    --input_fasta_file 'test-data/filter_alignment/alignment.fasta' \
    -o 'filter_alignment_default' \
    --allowed_gap_frac '0.999999' \
    --threshold '3.0'

filter_alignment.py \
    --input_fasta_file 'test-data/filter_alignment/alignment.fasta' \
    -o 'filter_alignment_without_mask_filter_and_outliers' \
    --suppress_lane_mask_filter \
    --allowed_gap_frac '0.999999' \
    --remove_outliers \
    --threshold '3.0'

filter_alignment.py \
    --input_fasta_file 'test-data/filter_alignment/alignment.fasta' \
    -o 'filter_alignment_entropy' \
    --allowed_gap_frac '0.999999' \
    --threshold '3.0' \
    --entropy_threshold '0.1'

# filter_fasta
filter_fasta.py \
    --input_fasta_fp 'test-data/filter_fasta/inseqs.fasta' \
    --output_fasta_fp 'filter_fasta_otu_map.fasta' \
    --otu_map 'test-data/filter_fasta/otu_map.txt'

filter_fasta.py \
    --input_fasta_fp 'test-data/filter_fasta/inseqs.fasta' \
    --output_fasta_fp 'filter_fasta_otu_map_negate.fasta' \
    --otu_map 'test-data/filter_fasta/otu_map.txt' \
    --negate

filter_fasta.py \
    --input_fasta_fp 'test-data/filter_fasta/inseqs.fasta' \
    --output_fasta_fp 'filter_fasta_seq_id.fasta' \
    --seq_id_fp 'test-data/filter_fasta/seqs_to_keep.txt'

filter_fasta.py \
    --input_fasta_fp 'test-data/filter_fasta/inseqs.fasta' \
    --output_fasta_fp 'filter_fasta_otu_table.fasta' \
    --biom_fp 'test-data/filter_fasta/otu_table.biom'

filter_fasta.py \
    --input_fasta_fp 'test-data/filter_fasta/inseqs.fasta' \
    --output_fasta_fp 'filter_fasta_subject_fasta.fasta' \
    --subject_fasta_fp 'test-data/filter_fasta/sl_inseqs.fasta'

filter_fasta.py \
    --input_fasta_fp 'test-data/filter_fasta/inseqs.fasta' \
    --output_fasta_fp 'filter_fasta_seq_id_prefix.fasta' \
    --seq_id_prefix 'S5'

filter_fasta.py \
    --input_fasta_fp 'test-data/filter_fasta/inseqs.fasta' \
    --output_fasta_fp 'filter_fasta_sample_id.fasta' \
    --sample_id_fp 'test-data/filter_fasta/map.txt'

# filter_otus_from_otu_table
filter_otus_from_otu_table.py \
    --input_fp 'test-data/filter_otus_from_otu_table/otu_table.biom' \
    --min_count '2' \
    --max_count '1000' \
    --min_samples '5' \
    --max_samples '350' \
    --output_fp 'test-data/filter_otus_from_otu_table/filtered_otu_table.biom'

filter_otus_from_otu_table.py \
    --input_fp 'test-data/filter_otus_from_otu_table/otu_table.biom' \
    --otu_ids_to_exclude_fp 'test-data/filter_otus_from_otu_table/chimeric_otus.txt' \
    --output_fp 'test-data/filter_otus_from_otu_table/chimera_filtered_otu_table.biom'

filter_otus_from_otu_table.py \
    --input_fp 'test-data/filter_otus_from_otu_table/otu_table.biom' \
    --otu_ids_to_exclude_fp 'test-data/filter_otus_from_otu_table/chimeric_otus.txt' \
    --negate_ids_to_exclude \
    --output_fp 'test-data/filter_otus_from_otu_table/chimera_picked_otu_table.biom'

# filter_samples_from_otu_table
filter_samples_from_otu_table.py \
    --input_fp 'test-data/filter_samples_from_otu_table/otu_table.biom' \
    --output_fp 'test-data/filter_samples_from_otu_table/tmp.biom' \
    --min_count '150'
biom convert \
    -i 'test-data/filter_samples_from_otu_table/tmp.biom' \
    -o 'test-data/filter_samples_from_otu_table/abundance_min.biom' \
    --to-json
rm 'test-data/filter_samples_from_otu_table/tmp.biom'

filter_samples_from_otu_table.py \
    --input_fp 'test-data/filter_samples_from_otu_table/otu_table.biom' \
    --output_fp 'test-data/filter_samples_from_otu_table/tmp.biom' \
    --min_count '0' \
    --max_count '149'
biom convert \
    -i 'test-data/filter_samples_from_otu_table/tmp.biom' \
    -o 'test-data/filter_samples_from_otu_table/abundance_max.biom' \
    --to-json
rm 'test-data/filter_samples_from_otu_table/tmp.biom'

filter_samples_from_otu_table.py \
    --input_fp 'test-data/filter_samples_from_otu_table/otu_table.biom' \
    --output_fp 'test-data/filter_samples_from_otu_table/tmp.biom' \
    --mapping_fp 'test-data/filter_samples_from_otu_table/map.txt' \
    --output_mapping_fp 'test-data/filter_samples_from_otu_table/metadata_positive.txt' \
    -s 'Treatment:Control'
biom convert \
    -i 'test-data/filter_samples_from_otu_table/tmp.biom' \
    -o 'test-data/filter_samples_from_otu_table/metadata_positive.biom' \
    --to-json
rm 'test-data/filter_samples_from_otu_table/tmp.biom'

filter_samples_from_otu_table.py \
    --input_fp 'test-data/filter_samples_from_otu_table/otu_table.biom' \
    --output_fp 'test-data/filter_samples_from_otu_table/tmp.biom' \
    --mapping_fp 'test-data/filter_samples_from_otu_table/map.txt' \
    -s 'Treatment:*,!Control'
biom convert \
    -i 'test-data/filter_samples_from_otu_table/tmp.biom' \
    -o 'test-data/filter_samples_from_otu_table/metadata_negative.biom' \
    --to-json
rm 'test-data/filter_samples_from_otu_table/tmp.biom'

filter_samples_from_otu_table.py \
    --input_fp 'test-data/filter_samples_from_otu_table/otu_table.biom' \
    --output_fp 'test-data/filter_samples_from_otu_table/tmp.biom' \
    --sample_id_fp 'test-data/filter_samples_from_otu_table/ids.txt'
biom convert \
    -i 'test-data/filter_samples_from_otu_table/tmp.biom' \
    -o 'test-data/filter_samples_from_otu_table/id_positive.biom' \
    --to-json
rm 'test-data/filter_samples_from_otu_table/tmp.biom'

filter_samples_from_otu_table.py \
    --input_fp 'test-data/filter_samples_from_otu_table/otu_table.biom' \
    --output_fp 'test-data/filter_samples_from_otu_table/tmp.biom' \
    --sample_id_fp 'test-data/filter_samples_from_otu_table/ids.txt' \
    --negate_sample_id_fp
biom convert \
    -i 'test-data/filter_samples_from_otu_table/tmp.biom' \
    -o 'test-data/filter_samples_from_otu_table/id_negative.biom' \
    --to-json
rm 'test-data/filter_samples_from_otu_table/tmp.biom'

# filter_taxa_from_otu_table
filter_taxa_from_otu_table.py \
    --input_otu_table_fp 'test-data/filter_taxa_from_otu_table/otu_table.biom' \
    --output_otu_table_fp 'test-data/filter_taxa_from_otu_table/tmp.biom' \
    --positive_taxa 'p__Bacteroidetes,p__Firmicutes' \
    --metadata_field 'taxonomy'
biom convert \
    -i 'test-data/filter_taxa_from_otu_table/tmp.biom' \
    -o 'test-data/filter_taxa_from_otu_table/positive_taxa.biom' \
    --to-json
rm 'test-data/filter_taxa_from_otu_table/tmp.biom'

filter_taxa_from_otu_table.py \
    --input_otu_table_fp 'test-data/filter_taxa_from_otu_table/otu_table.biom' \
    --output_otu_table_fp 'test-data/filter_taxa_from_otu_table/tmp.biom' \
    --negative_taxa 'p__Bacteroidetes,p__Firmicutes' \
    --metadata_field 'taxonomy'
biom convert \
    -i 'test-data/filter_taxa_from_otu_table/tmp.biom' \
    -o 'test-data/filter_taxa_from_otu_table/negative_taxa.biom' \
    --to-json
rm 'test-data/filter_taxa_from_otu_table/tmp.biom'

filter_taxa_from_otu_table.py \
    --input_otu_table_fp 'test-data/filter_taxa_from_otu_table/otu_table.biom' \
    --output_otu_table_fp 'test-data/filter_taxa_from_otu_table/tmp.biom' \
    --positive_taxa 'p__Firmicutes' \
    --negative_taxa 'c__Clostridia' \
    --metadata_field 'taxonomy'
biom convert \
    -i 'test-data/filter_taxa_from_otu_table/tmp.biom' \
    -o 'test-data/filter_taxa_from_otu_table/positive_negative_taxa.biom' \
    --to-json
rm 'test-data/filter_taxa_from_otu_table/tmp.biom'

# jackknifed_beta_diversity
jackknifed_beta_diversity.py \
    --otu_table_fp 'test-data/jackknifed_beta_diversity/otu_table.biom' \
    --mapping_fp 'test-data/jackknifed_beta_diversity/map.txt' \
    -o jackknifed_beta_diversity \
    --seqs_per_sample '10' \
    --tree_fp 'test-data/jackknifed_beta_diversity/rep_set.tre' \
    --master_tree 'consensus' \
    --parallel
rm -rf jackknifed_beta_diversity

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

# make_otu_heatmap
make_otu_heatmap.py \
    --otu_table_fp 'test-data/make_otu_heatmap/otu_table.biom' \
    --imagetype 'pdf' \
    --color_scheme "YlGn" \
    --width "5" \
    --height "5" \
    --dpi "200" \
    --obs_md_category "taxonomy" \
    --output_fp 'test-data/make_otu_heatmap/basic_heatmap.pdf'

make_otu_heatmap.py \
    --otu_table_fp 'test-data/make_otu_heatmap/otu_table.biom' \
    --imagetype 'png' \
    --color_scheme "YlGn" \
    --width "5" \
    --height "5" \
    --dpi "200" \
    --obs_md_category "taxonomy" \
    --output_fp 'test-data/make_otu_heatmap/basic_heatmap.png'

make_otu_heatmap.py \
    --otu_table_fp 'test-data/make_otu_heatmap/otu_table.biom' \
    --imagetype 'svg' \
    --color_scheme "YlGn" \
    --width "5" \
    --height "5" \
    --dpi "200" \
    --obs_md_category "taxonomy" \
    --output_fp 'test-data/make_otu_heatmap/basic_heatmap.svg'

make_otu_heatmap.py \
    --otu_table_fp 'test-data/make_otu_heatmap/otu_table.biom' \
    --map_fname 'test-data/make_otu_heatmap/mapping_file.txt' \
    --imagetype 'pdf' \
    --color_scheme "YlGn" \
    --width "5" \
    --height "5" \
    --dpi "200" \
    --obs_md_category "taxonomy" \
    --output_fp 'test-data/make_otu_heatmap/sample_sorted_heatmap.pdf'

make_otu_heatmap.py \
    --otu_table_fp 'test-data/make_otu_heatmap/otu_table.biom' \
    --map_fname 'test-data/make_otu_heatmap/mapping_file.txt' \
    --otu_tree 'test-data/make_otu_heatmap/rep_set.tre' \
    --imagetype 'pdf' \
    --color_scheme "YlGn" \
    --width "5" \
    --height "5" \
    --dpi "200" \
    --obs_md_category "taxonomy" \
    --output_fp 'test-data/make_otu_heatmap/sample_otu_sorted_heatmap.pdf'

make_otu_heatmap.py \
    --otu_table_fp 'test-data/make_otu_heatmap/otu_table.biom' \
    --map_fname 'test-data/make_otu_heatmap/mapping_file.txt' \
    --category "Treatment" \
    --imagetype 'pdf' \
    --color_scheme "YlGn" \
    --width "5" \
    --height "5" \
    --dpi "200" \
    --obs_md_category "taxonomy" \
    --output_fp 'test-data/make_otu_heatmap/treatment_sample_sorted_heatmap.pdf'

# make_phylogeny
make_phylogeny.py \
    --input_fp 'test-data/make_phylogeny/aligned.fasta' \
    --result_fp 'test-data/make_phylogeny/fasttree_tree_method_default.tre' \
    --tree_method 'fasttree' \
    --log_fp 'fasttree_tree_method_default.txt' \
    --root_method 'tree_method_default'

make_phylogeny.py \
    --input_fp 'test-data/make_phylogeny/aligned.fasta' \
    --result_fp 'raxml_v730.tre' \
    --tree_method 'raxml_v730' \
    --log_fp 'raxml_v730.txt' \
    --root_method 'tree_method_default'

make_phylogeny.py \
    --input_fp 'test-data/make_phylogeny/aligned.fasta' \
    --result_fp 'test-data/make_phylogeny/muscle.tre' \
    --tree_method 'muscle' \
    --log_fp 'muscle.txt' \
    --root_method 'tree_method_default'

make_phylogeny.py \
    --input_fp 'test-data/make_phylogeny/aligned.fasta' \
    --result_fp 'test-data/make_phylogeny/clustalw.tre' \
    --tree_method 'clustalw' \
    --log_fp 'clustalw.txt' \
    --root_method 'tree_method_default'

make_phylogeny.py \
    --input_fp 'test-data/make_phylogeny/aligned.fasta' \
    --result_fp 'clearcut.tre' \
    --tree_method 'clearcut' \
    --log_fp 'clearcut.txt' \
    --root_method 'tree_method_default'

make_phylogeny.py \
    --input_fp 'test-data/make_phylogeny/aligned.fasta' \
    --result_fp 'test-data/make_phylogeny/fasttree_midpoint.tre' \
    --tree_method 'fasttree' \
    --log_fp 'fasttree_midpoint.txt' \
    --root_method 'midpoint'

# multiple_join_paired_ends
multiple_join_paired_ends.py \
    --input_dir 'test-data/multiple_join_paired_ends/without_barcode/' \
    --output_dir 'test-data/multiple_join_paired_ends/output_without_barcode' \
    --read1_indicator 'forward_' \
    --read2_indicator 'reverse_' \
    --leading_text '' \
    --trailing_text ''

#multiple_join_paired_ends.py \
#    --input_dir 'test-data/multiple_join_paired_ends/without_barcode/' \
#    --output_dir 'multiple_join_paired_ends_without_barcode_parameter_files' \
#    --parameter_fp 'test-data/multiple_join_paired_ends/qiime_parameters.txt' \
#    --read1_indicator '_R1_' \
#    --read2_indicator '_R2_' \
#    --leading_text '' \
#    --trailing_text ''

multiple_join_paired_ends.py \
    --input_dir 'test-data/multiple_join_paired_ends/with_barcode/' \
    --output_dir 'test-data/multiple_join_paired_ends/output_with_barcode' \
    --read1_indicator 'forward_' \
    --read2_indicator 'reverse_' \
    --match_barcodes \
    --barcode_indicator 'barcode_' \
    --leading_text '' \
    --trailing_text ''

# multiple_split_libraries_fastq
multiple_split_libraries_fastq.py \
    --input_dir 'test-data/multiple_split_libraries_fastq/input' \
    --output_dir 'multiple_split_libraries_fastq' \
    --demultiplexing_method 'mapping_barcode_files' \
    --read_indicator 'reads_' \
    --barcode_indicator 'barcodes_' \
    --mapping_indicator 'mapping_' \
    --mapping_extensions 'txt' \
    --leading_text '' \
    --trailing_text '' \
    --sampleid_indicator '.'

multiple_split_libraries_fastq.py \
    --input_dir 'test-data/multiple_split_libraries_fastq/input' \
    --output_dir 'multiple_split_libraries_fastq_with_parameter_file' \
    --demultiplexing_method 'mapping_barcode_files' \
    --parameter_fp 'test-data/multiple_split_libraries_fastq/qiime_parameters.txt' \
    --read_indicator 'reads_' \
    --barcode_indicator 'barcodes_' \
    --mapping_indicator 'mapping_' \
    --mapping_extensions 'txt' \
    --leading_text '' \
    --trailing_text '' \
    --sampleid_indicator '.'

# pick_closed_reference_otus
pick_closed_reference_otus.py \
    --input_fp 'test-data/pick_closed_reference_otus/seqs.fna' \
    --output_dir 'pick_closed_reference_otus' \
    --reference_fp 'test-data/pick_closed_reference_otus/refseqs.fna' \
    --taxonomy_fp 'test-data/pick_closed_reference_otus/taxa.txt'
biom convert \
    -i 'pick_closed_reference_otus/otu_table.biom' \
    -o 'test-data/pick_closed_reference_otus/basic_otu_table.biom' \
    --to-json

pick_closed_reference_otus.py \
    --input_fp 'test-data/pick_closed_reference_otus/seqs.fna' \
    --output_dir 'pick_closed_reference_otus_sortmerna' \
    --reference_fp 'test-data/pick_closed_reference_otus/refseqs.fna' \
    --taxonomy_fp 'test-data/pick_closed_reference_otus/taxa.txt' \
    --parameter_fp 'test-data/pick_closed_reference_otus/sortmerna_params.txt'
biom convert \
    -i 'pick_closed_reference_otus_sortmerna/otu_table.biom' \
    -o 'test-data/pick_closed_reference_otus/sortmerna_otu_table.biom' \
    --to-json

pick_closed_reference_otus.py \
    --input_fp 'test-data/pick_closed_reference_otus/seqs.fna' \
    --output_dir 'pick_closed_reference_otus_assign_taxonomy' \
    --reference_fp 'test-data/pick_closed_reference_otus/refseqs.fna' \
    --assign_taxonomy
biom convert \
    -i 'pick_closed_reference_otus_assign_taxonomy/otu_table.biom' \
    -o 'test-data/pick_closed_reference_otus/assign_taxonomy_otu_table.biom' \
    --to-json

pick_closed_reference_otus.py \
    --input_fp 'test-data/pick_closed_reference_otus/seqs.fna' \
    --output_dir 'pick_closed_reference_otus_suppress_taxonomy_assignment' \
    --reference_fp 'test-data/pick_closed_reference_otus/refseqs.fna' \
    --suppress_taxonomy_assignment
biom convert \
    -i 'pick_closed_reference_otus_suppress_taxonomy_assignment/otu_table.biom' \
    -o 'test-data/pick_closed_reference_otus/suppress_taxonomy_assignment_otu_table.biom' \
    --to-json

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

# pick_otus
pick_otus.py \
    -i 'test-data/pick_otus/seqs.fna' \
    -o 'pick_otus_uclust' \
    --otu_picking_method 'uclust' \
    --similarity "0.97" \
    --denovo_otu_id_prefix "denovo" \
    --max_accepts "1" \
    --max_rejects "8" \
    --stepwords "8" \
    --word_length "8" \
    --non_chimeras_retention "union"

pick_otus.py \
    -i 'test-data/pick_otus/seqs.fna' \
    -o 'pick_otus_sortmerna' \
    --otu_picking_method "sortmerna" \
    --refseqs_fp "test-data/pick_otus/refseqs.fasta" \
    --sortmerna_e_value "1" \
    --sortmerna_coverage "0.97" \
    --sortmerna_tabular \
    --sortmerna_best_N_alignments "1" \
    --sortmerna_max_pos "10000" \
    --similarity "0.97" \
    --non_chimeras_retention "union"

#pick_otus.py \
#    -i 'test-data/pick_otus/seqs.fna' \
#    -o 'pick_otus_mothur' \
#    --otu_picking_method "mothur" \
#    --clustering_algorithm "furthest" \
#    --non_chimeras_retention "union"

pick_otus.py \
    -i 'test-data/pick_otus/seqs.fna' \
    -o 'pick_otus_trie' \
    --otu_picking_method "trie" \
    --non_chimeras_retention "union"

pick_otus.py \
    -i 'test-data/pick_otus/seqs.fna' \
    -o 'pick_otus_uclust_ref' \
    --otu_picking_method "uclust_ref" \
    --refseqs_fp "test-data/pick_otus/refseqs.fasta" \
    --similarity "0.97" \
    --max_accepts "1" \
    --max_rejects "8" \
    --stepwords "8" \
    --word_length "8" \
    --non_chimeras_retention "union"

#pick_otus.py \
#    -i 'test-data/pick_otus/seqs.fna' \
#    -o 'pick_otus_blast' \
#    --otu_picking_method "blast" \
#    --refseqs_fp "test-data/pick_otus/refseqs.fasta" \
#    --similarity "0.97" \
#    --max_e_value_blast "1e-10" \
#    --min_aligned_percent "0.5" \
#    --non_chimeras_retention "union"

# pick_otus.py \
#     -i 'test-data/pick_otus/seqs.fna' \
#     -o 'pick_otus_sumaclust' \
#     --otu_picking_method "sumaclust" \
#     --similarity "0.97" \
#     --sumaclust_l \
#     --denovo_otu_id_prefix "denovo" \
#     --non_chimeras_retention "union"

pick_otus.py \
    -i 'test-data/pick_otus/seqs.fna' \
    -o 'pick_otus_swarm' \
    --otu_picking_method "swarm" \
    --denovo_otu_id_prefix "denovo" \
    --swarm_resolution "1" \
    --non_chimeras_retention "union"

pick_otus.py \
    -i 'test-data/pick_otus/seqs.fna' \
    -o 'pick_otus_prefix_suffix' \
    --otu_picking_method "prefix_suffix" \
    --prefix_length "50" \
    --suffix_length "50" \
    --non_chimeras_retention "union"

pick_otus.py \
    -i 'test-data/pick_otus/seqs.fna' \
    -o 'pick_otus_cdhit' \
    --otu_picking_method "cdhit" \
    --similarity "0.97" \
    --non_chimeras_retention "union"

pick_otus.py \
    -i 'test-data/pick_otus/seqs.fna' \
    -o 'pick_otus_uclust_intersection' \
    --otu_picking_method "uclust" \
    --similarity "0.97" \
    --denovo_otu_id_prefix "denovo" \
    --max_accepts "1" \
    --max_rejects "8" \
    --stepwords "8" \
    --word_length "8" \
    --non_chimeras_retention "intersection"

# pick_rep_set
pick_rep_set.py \
    --input_file 'test-data/pick_rep_set/seqs_otus.txt' \
    --fasta_file 'test-data/pick_rep_set/seqs.fna' \
    --rep_set_picking_method 'first' \
    --sort_by 'otu' \
    --result_fp 'test-data/pick_rep_set/first_otu_fasta.fasta' \
    --log_fp 'test-data/pick_rep_set/first_otu_fasta.txt'

pick_rep_set.py \
    --input_file 'test-data/pick_rep_set/seqs_otus.txt' \
    --fasta_file 'test-data/pick_rep_set/seqs.fna' \
    --reference_seqs_fp 'test-data/pick_rep_set/refseqs.fasta' \
    --rep_set_picking_method 'first' \
    --sort_by 'otu' \
    --result_fp 'test-data/pick_rep_set/first_otu_fasta_ref.fasta' \
    --log_fp 'test-data/pick_rep_set/first_otu_fasta_ref.txt'

pick_rep_set.py \
    --input_file 'test-data/pick_rep_set/seqs_otus.txt' \
    --fasta_file 'test-data/pick_rep_set/seqs.fna' \
    --rep_set_picking_method 'longest' \
    --sort_by 'otu' \
    --result_fp 'test-data/pick_rep_set/longest_otu_fasta.fasta' \
    --log_fp 'test-data/pick_rep_set/longest_otu_fasta.txt'

pick_rep_set.py \
    --input_file 'test-data/pick_rep_set/seqs_otus.txt' \
    --fasta_file 'test-data/pick_rep_set/seqs.fna' \
    --rep_set_picking_method 'most_abundant' \
    --sort_by 'otu' \
    --result_fp 'test-data/pick_rep_set/most_abundant_otu_fasta.fasta' \
    --log_fp 'test-data/pick_rep_set/most_abundant_otu_fasta.txt'

pick_rep_set.py \
    --input_file 'test-data/pick_rep_set/seqs_otus.txt' \
    --fasta_file 'test-data/pick_rep_set/seqs.fna' \
    --rep_set_picking_method 'random' \
    --sort_by 'otu' \
    --result_fp 'test-data/pick_rep_set/random_otu_fasta.fasta' \
    --log_fp 'test-data/pick_rep_set/random_otu_fasta.txt'

pick_rep_set.py \
    --input_file 'test-data/pick_rep_set/seqs_otus.txt' \
    --fasta_file 'test-data/pick_rep_set/seqs.fna' \
    --rep_set_picking_method 'first' \
    --sort_by 'seq_id' \
    --result_fp 'test-data/pick_rep_set/first_seq_id_fasta.fasta' \
    --log_fp 'test-data/pick_rep_set/first_seq_id_fasta.txt'

# plot_taxa_summary
plot_taxa_summary.py \
   --counts_fname 'test-data/plot_taxa_summary/phylum.txt' \
   --dir_path 'test-data/plot_taxa_summary/phylum' \
   --labels 'phylum' \
   --num_categories '20' \
   --background_color 'white' \
   --dpi '80' \
   --x_width '12' \
   --y_height '12' \
   --bar_width '0.75' \
   --type_of_file 'png' \
   --chart_type 'area,bar,pie' \
   --resize_nth_label '0' \
   --label_type 'categorical'

plot_taxa_summary.py \
   --counts_fname 'test-data/plot_taxa_summary/phylum.txt,test-data/plot_taxa_summary/class.txt,test-data/plot_taxa_summary/genus.txt' \
   --dir_path 'test-data/plot_taxa_summary/phylum_class_genus' \
   --labels 'Phylum,Class,Genus' \
   --num_categories '20' \
   --background_color 'white' \
   --dpi '80' \
   --x_width '12' \
   --y_height '12' \
   --bar_width '0.75' \
   --type_of_file 'png' \
   --chart_type 'area,bar,pie' \
   --resize_nth_label '0' \
   --label_type 'categorical'

plot_taxa_summary.py \
   --counts_fname 'test-data/plot_taxa_summary/class.txt' \
   --dir_path 'test-data/plot_taxa_summary/class' \
   --labels 'Class' \
   --num_categories '10' \
   --background_color 'white' \
   --dpi '80' \
   --x_width '12' \
   --y_height '12' \
   --bar_width '0.75' \
   --chart_type 'pie' \
   --type_of_file 'svg' \
   --include_html_legend \
   --resize_nth_label '0' \
   --label_type 'categorical'

plot_taxa_summary.py \
   --counts_fname 'test-data/plot_taxa_summary/class.txt' \
   --dir_path 'test-data/plot_taxa_summary/class_colorby' \
   --labels 'Class' \
   --num_categories '20' \
   --colorby 'PC.636,PC.635' \
   --background_color 'white' \
   --dpi '80' \
   --x_width '12' \
   --y_height '12' \
   --bar_width '0.75' \
   --type_of_file 'pdf' \
   --chart_type 'area,bar,pie' \
   --resize_nth_label '0' \
   --label_type 'categorical'

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

# split_libraries_fastq
split_libraries_fastq.py \
    --sequence_read_fps 'test-data/split_libraries_fastq/lane1_read1.fastq.gz' \
    -o split_libraries_1 \
    --mapping_fps 'test-data/split_libraries_fastq/map.txt' \
    --barcode_read_fps 'test-data/split_libraries_fastq/lane1_barcode.fastq.gz' \
    --max_bad_run_length 3 \
    --min_per_read_length_fraction 0.75 \
    --sequence_max_n 0 \
    --start_seq_id 0 \
    --rev_comp_mapping_barcodes \
    --phred_quality_threshold 19 \
    --barcode_type 'golay_12' \
    --max_barcode_errors 1.5
rm -rf split_libraries_1

split_libraries_fastq.py \
    --sequence_read_fps 'test-data/split_libraries_fastq/lane1_read1.fastq.gz' \
    -o split_libraries_2 \
    --mapping_fps 'test-data/split_libraries_fastq/map.txt' \
    --barcode_read_fps 'test-data/split_libraries_fastq/lane1_barcode.fastq.gz' \
    --store_qual_scores \
    --store_demultiplexed_fastq \
    --max_bad_run_length 3 \
    --min_per_read_length_fraction 0.75 \
    --sequence_max_n 0 \
    --start_seq_id 0 \
    --rev_comp_mapping_barcodes \
    --phred_quality_threshold 19 \
    --barcode_type 'golay_12' \
    --max_barcode_errors 1.5
rm -rf split_libraries_2

split_libraries_fastq.py \
    --sequence_read_fps 'test-data/split_libraries_fastq/lane1_read1.fastq.gz,test-data/split_libraries_fastq/lane2_read1.fastq.gz' \
    -o split_libraries_3 \
    --mapping_fps 'test-data/split_libraries_fastq/map.txt,test-data/split_libraries_fastq/map.txt' \
    --barcode_read_fps 'test-data/split_libraries_fastq/lane1_barcode.fastq.gz,test-data/split_libraries_fastq/lane2_barcode.fastq.gz' \
    --max_bad_run_length 3 \
    --min_per_read_length_fraction 0.75 \
    --sequence_max_n 0 \
    --start_seq_id 0 \
    --rev_comp_mapping_barcodes \
    --phred_quality_threshold 19 \
    --barcode_type 'golay_12' \
    --max_barcode_errors 1.5
rm -rf split_libraries_3

split_libraries_fastq.py \
    --sequence_read_fps 'test-data/split_libraries_fastq/lane1_read1.fastq.gz' \
    -o split_libraries_4 \
    --sample_ids 'my.sample.1' \
    --max_bad_run_length 3 \
    --min_per_read_length_fraction 0.75 \
    --sequence_max_n 0 \
    --start_seq_id 0 \
    --rev_comp_mapping_barcodes \
    --phred_quality_threshold 19 \
    --barcode_type 'not-barcoded' \
    --max_barcode_errors 1.5
rm -rf split_libraries_4

split_libraries_fastq.py \
    --sequence_read_fps 'test-data/split_libraries_fastq/lane1_read1.fastq.gz,test-data/split_libraries_fastq/lane2_read1.fastq.gz' \
    -o split_libraries_5 \
    --sample_ids 'my.sample.1,my.sample.2' \
    --max_bad_run_length 3 \
    --min_per_read_length_fraction 0.75 \
    --sequence_max_n 0 \
    --start_seq_id 0 \
    --rev_comp_mapping_barcodes \
    --phred_quality_threshold 19 \
    --barcode_type 'not-barcoded' \
    --max_barcode_errors 1.5
rm -rf split_libraries_5

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

# summarize_taxa_through_plots
summarize_taxa_through_plots.py \
    --otu_table_fp 'test-data/summarize_taxa_through_plots/otu_table.biom' \
    --output_dir summarize_taxa_through_plots_mapping \
    --mapping_fp 'test-data/summarize_taxa_through_plots/Fasting_Map.txt'
biom convert \
    -i 'summarize_taxa_through_plots_mapping/otu_table_L2.biom' \
    -o 'summarize_taxa_through_plots_mapping/otu_table_L2_json.biom' \
    --to-json
biom convert \
    -i 'summarize_taxa_through_plots_mapping/otu_table_L3.biom' \
    -o 'summarize_taxa_through_plots_mapping/otu_table_L3_json.biom' \
    --to-json
biom convert \
    -i 'summarize_taxa_through_plots_mapping/otu_table_L4.biom' \
    -o 'summarize_taxa_through_plots_mapping/otu_table_L4_json.biom' \
    --to-json
biom convert \
    -i 'summarize_taxa_through_plots_mapping/otu_table_L5.biom' \
    -o 'summarize_taxa_through_plots_mapping/otu_table_L5_json.biom' \
    --to-json
biom convert \
    -i 'summarize_taxa_through_plots_mapping/otu_table_L6.biom' \
    -o 'summarize_taxa_through_plots_mapping/otu_table_L6_json.biom' \
    --to-json
cp summarize_taxa_through_plots_mapping/*.txt test-data/summarize_taxa_through_plots/mapping/
cp summarize_taxa_through_plots_mapping/*_json.biom test-data/summarize_taxa_through_plots/mapping/
cp summarize_taxa_through_plots_mapping/taxa_summary_plots/area_charts.html 'test-data/summarize_taxa_through_plots/mapping/area_charts.html'
cp summarize_taxa_through_plots_mapping/taxa_summary_plots/bar_charts.html 'test-data/summarize_taxa_through_plots/mapping/bar_charts.html'

summarize_taxa_through_plots.py \
    --otu_table_fp 'test-data/summarize_taxa_through_plots/otu_table.biom' \
    --output_dir summarize_taxa_through_plots_mapping_categories \
    --mapping_fp 'test-data/summarize_taxa_through_plots/Fasting_Map.txt' \
    --mapping_category 'Treatment'
biom convert \
    -i 'summarize_taxa_through_plots_mapping_categories/Treatment_otu_table_L2.biom' \
    -o 'summarize_taxa_through_plots_mapping_categories/Treatment_otu_table_L2_json.biom' \
    --to-json
biom convert \
    -i 'summarize_taxa_through_plots_mapping_categories/Treatment_otu_table_L3.biom' \
    -o 'summarize_taxa_through_plots_mapping_categories/Treatment_otu_table_L3_json.biom' \
    --to-json
biom convert \
    -i 'summarize_taxa_through_plots_mapping_categories/Treatment_otu_table_L4.biom' \
    -o 'summarize_taxa_through_plots_mapping_categories/Treatment_otu_table_L4_json.biom' \
    --to-json
biom convert \
    -i 'summarize_taxa_through_plots_mapping_categories/Treatment_otu_table_L5.biom' \
    -o 'summarize_taxa_through_plots_mapping_categories/Treatment_otu_table_L5_json.biom' \
    --to-json
biom convert \
    -i 'summarize_taxa_through_plots_mapping_categories/Treatment_otu_table_L6.biom' \
    -o 'summarize_taxa_through_plots_mapping_categories/Treatment_otu_table_L6_json.biom' \
    --to-json
cp summarize_taxa_through_plots_mapping_categories/*.txt test-data/summarize_taxa_through_plots/mapping_categories/
cp summarize_taxa_through_plots_mapping_categories/*_json.biom test-data/summarize_taxa_through_plots/mapping_categories/
cp summarize_taxa_through_plots_mapping_categories/taxa_summary_plots/area_charts.html 'test-data/summarize_taxa_through_plots/mapping_categories/area_charts.html'
cp summarize_taxa_through_plots_mapping_categories/taxa_summary_plots/bar_charts.html 'test-data/summarize_taxa_through_plots/mapping_categories/bar_charts.html'

summarize_taxa_through_plots.py \
    --otu_table_fp 'test-data/summarize_taxa_through_plots/otu_table.biom' \
    --output_dir summarize_taxa_through_plots_mapping_sort \
    --mapping_fp 'test-data/summarize_taxa_through_plots/Fasting_Map.txt' \
    --sort
biom convert \
    -i 'summarize_taxa_through_plots_mapping_sort/otu_table_sorted_L2.biom' \
    -o 'summarize_taxa_through_plots_mapping_sort/otu_table_sorted_L2_json.biom' \
    --to-json
biom convert \
    -i 'summarize_taxa_through_plots_mapping_sort/otu_table_sorted_L3.biom' \
    -o 'summarize_taxa_through_plots_mapping_sort/otu_table_sorted_L3_json.biom' \
    --to-json
biom convert \
    -i 'summarize_taxa_through_plots_mapping_sort/otu_table_sorted_L4.biom' \
    -o 'summarize_taxa_through_plots_mapping_sort/otu_table_sorted_L4_json.biom' \
    --to-json
biom convert \
    -i 'summarize_taxa_through_plots_mapping_sort/otu_table_sorted_L5.biom' \
    -o 'summarize_taxa_through_plots_mapping_sort/otu_table_sorted_L5_json.biom' \
    --to-json
biom convert \
    -i 'summarize_taxa_through_plots_mapping_sort/otu_table_sorted_L6.biom' \
    -o 'summarize_taxa_through_plots_mapping_sort/otu_table_sorted_L6_json.biom' \
    --to-json
cp summarize_taxa_through_plots_mapping_sort/*.txt test-data/summarize_taxa_through_plots/mapping_sort/
cp summarize_taxa_through_plots_mapping_sort/*_json.biom test-data/summarize_taxa_through_plots/mapping_sort/
cp summarize_taxa_through_plots_mapping_sort/taxa_summary_plots/area_charts.html 'test-data/summarize_taxa_through_plots/mapping_sort/area_charts.html'
cp summarize_taxa_through_plots_mapping_sort/taxa_summary_plots/bar_charts.html 'test-data/summarize_taxa_through_plots/mapping_sort/bar_charts.html'

summarize_taxa_through_plots.py \
    --otu_table_fp 'test-data/summarize_taxa_through_plots/otu_table.biom' \
    --output_dir summarize_taxa_through_plots_without_mapping
biom convert \
    -i 'summarize_taxa_through_plots_without_mapping/otu_table_L2.biom' \
    -o 'summarize_taxa_through_plots_without_mapping/otu_table_L2_json.biom' \
    --to-json
biom convert \
    -i 'summarize_taxa_through_plots_without_mapping/otu_table_L3.biom' \
    -o 'summarize_taxa_through_plots_without_mapping/otu_table_L3_json.biom' \
    --to-json
biom convert \
    -i 'summarize_taxa_through_plots_without_mapping/otu_table_L4.biom' \
    -o 'summarize_taxa_through_plots_without_mapping/otu_table_L4_json.biom' \
    --to-json
biom convert \
    -i 'summarize_taxa_through_plots_without_mapping/otu_table_L5.biom' \
    -o 'summarize_taxa_through_plots_without_mapping/otu_table_L5_json.biom' \
    --to-json
biom convert \
    -i 'summarize_taxa_through_plots_without_mapping/otu_table_L6.biom' \
    -o 'summarize_taxa_through_plots_without_mapping/otu_table_L6_json.biom' \
    --to-json
cp summarize_taxa_through_plots_without_mapping/*.txt test-data/summarize_taxa_through_plots/without_mapping/
cp summarize_taxa_through_plots_without_mapping/*_json.biom test-data/summarize_taxa_through_plots/without_mapping/
cp summarize_taxa_through_plots_without_mapping/taxa_summary_plots/area_charts.html 'test-data/summarize_taxa_through_plots/without_mapping/area_charts.html'
cp summarize_taxa_through_plots_without_mapping/taxa_summary_plots/bar_charts.html 'test-data/summarize_taxa_through_plots/without_mapping/bar_charts.html'

# upgma_cluster
upgma_cluster.py \
    --input_path 'test-data/upgma_cluster/' \
    --output_path 'test-data/upgma_cluster/'

# validate_mapping_file
validate_mapping_file.py \
    -m 'test-data/validate_mapping_file/map.tsv' \
    -o validate_mapping_file_output \
    -c '_'
cp validate_mapping_file_output/*.html 'test-data/validate_mapping_file/map.tsv.html'
cp validate_mapping_file_output/*.log 'test-data/validate_mapping_file/map.tsv.log'
cp validate_mapping_file_output/*corrected.txt 'test-data/validate_mapping_file/map.tsv_corrected.txt'
rm -rf validate_mapping_file_output
