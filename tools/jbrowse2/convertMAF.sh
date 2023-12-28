#!/usr/bin/env bash
# https://github.com/cmdcolin/jbrowse-plugin-mafviewer/blob/master/bin/convert.sh
#  MAF file must contain the species name and chromosome name
#  e.g. hg38.chr1 in the sequence identifiers.
$3/maf2bed.pl $2 < $1 > `basename $1`.bed
sort -k1,1 -k2,2n `basename $1`.bed > `basename $1`.sorted.bed
bgzip `basename $1`.sorted.bed
tabix -p bed `basename $1`.sorted.bed.gz
