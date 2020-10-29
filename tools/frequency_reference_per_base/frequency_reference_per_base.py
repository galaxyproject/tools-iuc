#!/usr/bin/env python3

# Takes in genome annotation file and a bam alignment to that file
# Generates x4 bigWig files:
#   1. Frequency of aligned matches to reference
#   2. Frequency of aligned mismatches to reference
#   3. Frequency of aligned insertions to reference
#   4. Frequency of aligned deletions to reference
#
# Usage statement:
# python getPercentReferencePerBase.py in.bam ref.fasta prefix
#
# 10/28/2020 - Nathan P. Roach, natproach@gmail.com
import argparse

import numpy as np
import pyBigWig
import pysam
import pysamstats


parser = argparse.ArgumentParser(
    description='Generate bigWigs of frequencies'
    ' of matches, mismatches, insertions,'
    ' and deletions at each position in a reference genome')

parser.add_argument(
    '--bam',
    help='Alignment BAM to use to generate pileup')

parser.add_argument(
    '--reference',
    help='Reference genome FASTA to use to generate pileup')

parser.add_argument('--prefix', help="Prefix for output files")

parser.add_argument(
    '--depth', type=int,
    help='Maximum depth to consider in pileup step')

args = parser.parse_args()

bam_infilepath = args.bam
bam_in = pysam.AlignmentFile(bam_infilepath, "rb")

ref_infilepath = args.reference
pileup_stats = pysamstats.load_pileup("variation",
                                      bam_in,
                                      fafile=ref_infilepath,
                                      pad=True,
                                      max_depth=args.depth)

prefix = args.prefix

chroms = pileup_stats.chrom
starts = pileup_stats.pos
ends = np.add(pileup_stats.pos, 1)
match_percent = np.divide(pileup_stats.matches, pileup_stats.reads_all)
mismatch_percent = np.divide(pileup_stats.mismatches, pileup_stats.reads_all)
insertion_percent = np.divide(pileup_stats.insertions, pileup_stats.reads_all)
deletion_percent = np.divide(pileup_stats.deletions, pileup_stats.reads_all)

# Calculate header for bw file
bw_header = []
for i in range(bam_in.nreferences):
    chr_name = bam_in.get_reference_name(i)
    bw_header.append((chr_name, bam_in.get_reference_length(chr_name)))

output_suffixes = ["match", "mismatch", "insertion", "deletion"]
output_data = [match_percent,
               mismatch_percent,
               insertion_percent,
               deletion_percent]
for i in range(len(output_suffixes)):
    bw = pyBigWig.open(prefix + "_" + output_suffixes[i] + ".bw", "w")
    bw.addHeader(bw_header)
    bw.addEntries(chroms.astype("str"),
                  starts,
                  ends=ends,
                  values=output_data[i])
    bw.close()
