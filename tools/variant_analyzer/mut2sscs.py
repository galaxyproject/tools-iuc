#!/usr/bin/env python

"""mut2sscs.py

Author -- Gundula Povysil
Contact -- povysil@bioinf.jku.at

Takes a tabular file with mutations from DCS and a BAM file of SSCS as input
and extracts all tags of reads that carry the mutation.
Calculates statistics about number of ab/ba/duplex per mutation.

=======  ==========  =================  ================================
Version  Date        Author             Description
0.2.1    2019-10-27  Gundula Povysil    -
=======  ==========  =================  ================================

USAGE: python mut2sscs.py DCS_Mutations.tabular SSCS.bam SSCS_counts.json

"""

from __future__ import division

import argparse
import json
import os
import sys

import pysam
from cyvcf2 import VCF


def make_argparser():
    parser = argparse.ArgumentParser(description='Takes a vcf file with mutations and a BAM file as input and prints all tags of reads that carry the mutation to a user specified output file.')
    parser.add_argument('--mutFile',
                        help='VCR file with DCS mutations.')
    parser.add_argument('--bamFile',
                        help='BAM file with aligned SSCS reads.')
    parser.add_argument('--outputJson',
                        help='Output JSON file to store SSCS counts.')
    return parser


def mut2sscs(argv):
    parser = make_argparser()
    args = parser.parse_args(argv[1:])

    file1 = args.mutFile
    file2 = args.bamFile
    sscs_counts_json = args.outputJson

    if os.path.isfile(file1) is False:
        sys.exit("Error: Could not find '{}'".format(file1))

    if os.path.isfile(file2) is False:
        sys.exit("Error: Could not find '{}'".format(file2))

    # read SSCS bam file
#    pysam.index(file2)
    bam = pysam.AlignmentFile(file2, "rb")

    # get tags
    mut_pos_dict = {}
    ref_pos_dict = {}

    for variant in VCF(file1):
        chrom = variant.CHROM
        stop_pos = variant.start
        chrom_stop_pos = str(chrom) + "#" + str(stop_pos)
        ref = variant.REF
        alt = variant.ALT[0]
#        nc = variant.format('NC')
        # ad = variant.format('AD')

        if len(ref) == len(alt):

            for pileupcolumn in bam.pileup(chrom, stop_pos - 1, stop_pos + 1, max_depth=1000000000):
                if pileupcolumn.reference_pos == stop_pos:
                    count_alt = 0
                    count_ref = 0
                    count_indel = 0
                    print("unfiltered reads=", pileupcolumn.n, "filtered reads=", len(pileupcolumn.pileups),
                          "difference= ", len(pileupcolumn.pileups) - pileupcolumn.n)
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            tag = pileupread.alignment.query_name
                            abba = tag[-2:]
                            # query position is None if is_del or is_refskip is set.
                            if pileupread.alignment.query_sequence[pileupread.query_position] == alt:
                                count_alt += 1
                                if chrom_stop_pos in mut_pos_dict:
                                    if abba in mut_pos_dict[chrom_stop_pos]:
                                        mut_pos_dict[chrom_stop_pos][abba] += 1
                                    else:
                                        mut_pos_dict[chrom_stop_pos][abba] = 1
                                else:
                                    mut_pos_dict[chrom_stop_pos] = {}
                                    mut_pos_dict[chrom_stop_pos][abba] = 1
                                if chrom_stop_pos not in ref_pos_dict:
                                    ref_pos_dict[chrom_stop_pos] = {}
                                    ref_pos_dict[chrom_stop_pos][abba] = 0

                            elif pileupread.alignment.query_sequence[pileupread.query_position] == ref:
                                count_ref += 1
                                if chrom_stop_pos in ref_pos_dict:
                                    if abba in ref_pos_dict[chrom_stop_pos]:
                                        ref_pos_dict[chrom_stop_pos][abba] += 1
                                    else:
                                        ref_pos_dict[chrom_stop_pos][abba] = 1
                                else:
                                    ref_pos_dict[chrom_stop_pos] = {}
                                    ref_pos_dict[chrom_stop_pos][abba] = 1
                            else:
                                count_indel += 1

                    print("coverage at pos %s = %s, ref = %s, alt = %s, indel = %s,\n" %
                          (pileupcolumn.pos, count_ref + count_alt, count_ref, count_alt, count_indel))

            # if mutation is in DCS file but not in SSCS, then set counts to NA
            if chrom_stop_pos not in mut_pos_dict.keys():
                mut_pos_dict[chrom_stop_pos] = {}
                mut_pos_dict[chrom_stop_pos]["ab"] = 0
                mut_pos_dict[chrom_stop_pos]["ba"] = 0
                ref_pos_dict[chrom_stop_pos] = {}
                ref_pos_dict[chrom_stop_pos]["ab"] = 0
                ref_pos_dict[chrom_stop_pos]["ba"] = 0
        else:
            print("indels are currently not evaluated")
    bam.close()

    # save counts
    with open(sscs_counts_json, "w") as f:
        json.dump((mut_pos_dict, ref_pos_dict), f)


if __name__ == '__main__':
    sys.exit(mut2sscs(sys.argv))
