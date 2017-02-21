#!/usr/bin/env bash

# make_otu_table
make_otu_table.py \
    --otu_map_fp 'test-data/make_otu_table/otu_map.txt' \
    --taxonomy 'test-data/make_otu_table/tax_assignments.txt' \
    --exclude_otus_fp 'test-data/make_otu_table/chimeric_seqs.txt' \
    --mapping_fp 'test-data/make_otu_table/mapping_file.txt' \
    --output_biom_fp 'test-data/make_otu_table/OTU_table_chimeric.biom'

make_otu_table.py \
    --otu_map_fp 'test-data/make_otu_table/otu_map.txt' \
    --taxonomy 'test-data/make_otu_table/tax_assignments.txt' \
    --exclude_otus_fp 'test-data/make_otu_table/pynast_failures.fna' \
    --mapping_fp 'test-data/make_otu_table/mapping_file.txt' \
    --output_biom_fp 'test-data/make_otu_table/OTU_table_pynast.biom'
