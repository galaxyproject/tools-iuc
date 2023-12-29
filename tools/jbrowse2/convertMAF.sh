#!/usr/bin/env bash
# https://github.com/cmdcolin/jbrowse-plugin-mafviewer/blob/master/bin/convert.sh
#  MAF file must contain the species name and chromosome name
#  e.g. hg38.chr1 in the sequence identifiers.
perl $3/maf2bed.pl $2 < $1 > $1.bed
sort -k1,1 -k2,2n $1.bed > $1.sorted.bed
bgzip -c $1.sorted.bed > $1.sorted.bed.gz
tabix -p bed $1.sorted.bed.gz
