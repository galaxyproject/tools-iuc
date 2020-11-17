#!/usr/bin/env python

"""mut2read.py

Author -- Gundula Povysil
Contact -- povysil@bioinf.jku.at

Takes a tabular file with mutations and a BAM file as input and prints
all tags of reads that carry the mutation to a user specified output file.
Creates fastq file of reads of tags with mutation.

=======  ==========  =================  ================================
Version  Date        Author             Description
2.0.0    2020-10-30  Gundula Povysil    -
=======  ==========  =================  ================================

USAGE: python mut2read.py DCS_Mutations.tabular DCS.bam Aligned_Families.tabular Interesting_Reads.fastq tag_count_dict.json
"""

import argparse
import json
import os
import sys

import numpy as np
import pysam
from cyvcf2 import VCF


def make_argparser():
    parser = argparse.ArgumentParser(description='Takes a vcf file with mutations and a BAM file as input and prints all tags of reads that carry the mutation to a user specified output file and creates a fastq file of reads of tags with mutation.')
    parser.add_argument('--mutFile',
                        help='VCF file with DCS mutations.')
    parser.add_argument('--bamFile',
                        help='BAM file with aligned DCS reads.')
    parser.add_argument('--familiesFile',
                        help='TABULAR file with aligned families.')
    parser.add_argument('--outputFastq',
                        help='Output FASTQ file of reads with mutations.')
    parser.add_argument('--outputJson',
                        help='Output JSON file to store collected data.')
    return parser


def mut2read(argv):
    parser = make_argparser()
    args = parser.parse_args(argv[1:])

    file1 = args.mutFile
    file2 = args.bamFile
    file3 = args.familiesFile
    outfile = args.outputFastq
    json_file = args.outputJson

    if os.path.isfile(file1) is False:
        sys.exit("Error: Could not find '{}'".format(file1))

    if os.path.isfile(file2) is False:
        sys.exit("Error: Could not find '{}'".format(file2))

    if os.path.isfile(file3) is False:
        sys.exit("Error: Could not find '{}'".format(file3))

    # read dcs bam file
    bam = pysam.AlignmentFile(file2, "rb")

    # get tags
    tag_dict = {}
    cvrg_dict = {}

    for variant in VCF(file1):
        chrom = variant.CHROM
        stop_pos = variant.start
        chrom_stop_pos = str(chrom) + "#" + str(stop_pos)
        ref = variant.REF
        alt = variant.ALT[0]
        dcs_len = []
        if len(ref) == len(alt):
            for pileupcolumn in bam.pileup(chrom, stop_pos - 1, stop_pos + 1, max_depth=100000000):
                if pileupcolumn.reference_pos == stop_pos:
                    count_alt = 0
                    count_ref = 0
                    count_indel = 0
                    count_n = 0
                    count_other = 0
                    count_lowq = 0
                    print("unfiltered reads=", pileupcolumn.n, "filtered reads=", len(pileupcolumn.pileups),
                          "difference= ", len(pileupcolumn.pileups) - pileupcolumn.n)
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            # query position is None if is_del or is_refskip is set.
                            nuc = pileupread.alignment.query_sequence[pileupread.query_position]
                            dcs_len.append(len(pileupread.alignment.query_sequence))
                            if nuc == alt:
                                count_alt += 1
                                tag = pileupread.alignment.query_name
                                if tag in tag_dict:
                                    tag_dict[tag][chrom_stop_pos] = alt
                                else:
                                    tag_dict[tag] = {}
                                    tag_dict[tag][chrom_stop_pos] = alt
                            elif nuc == ref:
                                count_ref += 1
                            elif nuc == "N":
                                count_n += 1
                            elif nuc == "lowQ":
                                count_lowq += 1
                            else:
                                count_other += 1
                        else:
                            count_indel += 1
                    dcs_median = np.median(np.array(dcs_len))
                    cvrg_dict[chrom_stop_pos] = (count_ref, count_alt, dcs_median)

                    print("coverage at pos %s = %s, ref = %s, alt = %s, other bases = %s, N = %s, indel = %s, low quality = %s, median length of DCS = %s\n" %
                          (pileupcolumn.pos, count_ref + count_alt, count_ref, count_alt, count_other, count_n,
                           count_indel, count_lowq, dcs_median))
        else:
            print("indels are currently not evaluated")
    bam.close()

    with open(json_file, "w") as f:
        json.dump((tag_dict, cvrg_dict), f)

    # create fastq from aligned reads
    with open(outfile, 'w') as out:
        with open(file3, 'r') as families:
            for line in families:
                line = line.rstrip('\n')
                splits = line.split('\t')
                tag = splits[0]

                if tag in tag_dict:
                    str1 = splits[4]
                    curr_seq = str1.replace("-", "")
                    str2 = splits[5]
                    curr_qual = str2.replace(" ", "")

                    out.write("@" + splits[0] + "." + splits[1] + "." + splits[2] + "\n")
                    out.write(curr_seq + "\n")
                    out.write("+" + "\n")
                    out.write(curr_qual + "\n")


if __name__ == '__main__':
    sys.exit(mut2read(sys.argv))
