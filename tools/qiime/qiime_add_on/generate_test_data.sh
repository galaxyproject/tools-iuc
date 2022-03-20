#!/usr/bin/env bash

# make_otu_table
make_otu_table.py \
    --otu_map_fp 'test-data/make_otu_table/otu_map.txt' \
    --taxonomy 'test-data/make_otu_table/tax_assignments.txt' \
    --exclude_otus_fp 'test-data/make_otu_table/chimeric_seqs.txt' \
    --mapping_fp 'test-data/make_otu_table/mapping_file.txt' \
    --output_biom_fp 'test-data/make_otu_table/OTU_table_chimeric.biom'
biom convert \
    -i 'test-data/make_otu_table/OTU_table_chimeric.biom' \
    -o 'test-data/make_otu_table/OTU_table_chimeric.biom' \
    --to-json

make_otu_table.py \
    --otu_map_fp 'test-data/make_otu_table/otu_map.txt' \
    --taxonomy 'test-data/make_otu_table/tax_assignments.txt' \
    --exclude_otus_fp 'test-data/make_otu_table/pynast_failures.fna' \
    --mapping_fp 'test-data/make_otu_table/mapping_file.txt' \
    --output_biom_fp 'test-data/make_otu_table/OTU_table_pynast.biom'
biom convert \
    -i 'test-data/make_otu_table/OTU_table_pynast.biom' \
    -o 'test-data/make_otu_table/OTU_table_pynast.biom' \
    --to-json

# collapse_samples
collapse_samples.py \
    --input_biom_fp 'test-data/collapse_samples/table.biom' \
    --mapping_fp 'test-data/collapse_samples/map.txt' \
    --collapse_mode 'sum' \
    --collapse_fields 'SampleType' \
    --output_biom_fp 'test-data/collapse_samples/collapsed_sum_SampleType_table.biom' \
    --output_mapping_fp 'test-data/collapse_samples/collapsed_sum_SampleType_map.txt'
biom convert \
    -i 'test-data/collapse_samples/collapsed_sum_SampleType_table.biom' \
    -o 'test-data/collapse_samples/collapsed_sum_SampleType_table.biom' \
    --to-json
    
collapse_samples.py \
    --input_biom_fp 'test-data/collapse_samples/table.biom' \
    --mapping_fp 'test-data/collapse_samples/map.txt' \
    --collapse_mode 'first' \
    --collapse_fields 'subject','year' \
    --normalize \
    --output_biom_fp 'test-data/collapse_samples/collapsed_first_2fields_table.biom' \
    --output_mapping_fp 'test-data/collapse_samples/collapsed_first_2fields_map.txt'
biom convert \
    -i 'test-data/collapse_samples/collapsed_first_2fields_table.biom' \
    -o 'test-data/collapse_samples/collapsed_first_2fields_table.biom' \
    --to-json