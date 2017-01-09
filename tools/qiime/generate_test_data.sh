#!/usr/bin/env bash

# validate_mapping_file
validate_mapping_file.py \
    -m 'test-data/map.tsv' \
    -o validate_mapping_file_output \
    -c '_'
cp validate_mapping_file_output/*.html 'test-data/map.tsv.html'
cp validate_mapping_file_output/*.log 'test-data/map.tsv.log'
cp validate_mapping_file_output/*corrected.txt 'test-data/map.tsv_corrected.txt'
rm -rf validate_mapping_file_output

