#!/usr/bin/env bash
# https://github.com/cmdcolin/jbrowse-plugin-mafviewer/blob/master/bin/convert.sh
# maf2bed modified to work right as a python script by ross lazarus in desperation
#  MAF file must contain the species name and chromosome name
#  e.g. hg38.chr1 in the sequence identifiers with hg38 passed in as $2
python $3/maf2bed.py $2 < $1 | sort -k1,1 -k2,2n > $4.sorted.bed
bgzip -c $4.sorted.bed > $4.sorted.bed.gz
tabix -p bed $4.sorted.bed.gz
