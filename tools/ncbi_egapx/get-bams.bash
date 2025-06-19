#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail

##
## Collect output bam files from STAR
##

bam_temp_dir=$(mktemp --directory ./egapx_out/bam.XXXXXXXXXX) || {
    echo "ERROR: failed making temporary directory"
    exit 1
}

bam_files_file="$bam_temp_dir/bam-files.txt"
touch "$bam_files_file" || {
    echo "ERROR: failed creating temporary file"
    exit 1
}
trap 'rm -f -- "$bam_files_file"' EXIT

find . -type f -name \*Aligned.out.Sorted.bam -print > "$bam_files_file" || {
    echo "ERROR: failed finding bam files"
    exit 1
}

idx=0
while IFS="" read pathname; do
    filename=$(basename "$pathname" ".bam")

    dest_path_prefix="$bam_temp_dir/$filename"
    if [ -e "${dest_path_prefix}.bam" ]; then
        dest_path_prefix="${dest_path_prefix}.${idx}"
        idx=$((idx + 1))
    fi
    dest_pathname="${dest_path_prefix}.bam"

    cp "$pathname" "$dest_pathname" || {
        echo "ERROR: failed: cp \"$pathname\" \"$dest_pathname\""
        exit 1
    }

done < "$bam_files_file"

exit 0
