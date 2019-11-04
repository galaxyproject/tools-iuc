#!/usr/bin/env python

"""mut2read.py

Author -- Gundula Povysil
Contact -- povysil@bioinf.jku.at

Takes a tabular file with mutations and a BAM file as input and prints
all tags of reads that carry the mutation to a user specified output file.
Creates fastq file of reads of tags with mutation.

=======  ==========  =================  ================================
Version  Date        Author             Description
0.2.1    2019-10-27  Gundula Povysil    -
=======  ==========  =================  ================================

USAGE: python mut2read.py DCS_Mutations.tabular DCS.bam Aligned_Families.tabular Interesting_Reads.fastq
                          tag_count_dict.json
"""

import argparse
import json
import os
import sys

import numpy as np
import pysam


def make_argparser():
    parser = argparse.ArgumentParser(description='Takes a tabular file with mutations and a BAM file as input and prints all tags of reads that carry the mutation to a user specified output file and creates a fastq file of reads of tags with mutation.')
    parser.add_argument('--mutFile',
                        help='TABULAR file with DCS mutations.')
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
        print("Error: Could not find '{}'".format(file1))
        exit(0)

    if os.path.isfile(file2) is False:
        print("Error: Could not find '{}'".format(file2))
        exit(1)

    if os.path.isfile(file3) is False:
        print("Error: Could not find '{}'".format(file3))
        exit(2)

    # read mut file
    with open(file1, 'r') as mut:
        mut_array = np.genfromtxt(mut, skip_header=1, delimiter='\t', comments='#', dtype='string')

    # read dcs bam file
    pysam.index(file2)
    bam = pysam.AlignmentFile(file2, "rb")

    # get tags
    mut_tags = []
    tag_dict = {}
    cvrg_dict = {}

    if len(mut_array) == 13:
        mut_array = mut_array.reshape((1, len(mut_array)))

    for m in range(len(mut_array[:, 0])):
        print(str(m + 1) + " of " + str(len(mut_array[:, 0])))
        chrom = mut_array[m, 1]
        stop_pos = mut_array[m, 2].astype(int)
        chrom_stop_pos = str(chrom) + "#" + str(stop_pos)
        ref = mut_array[m, 9]
        alt = mut_array[m, 10]

        dcs_len = []

        for pileupcolumn in bam.pileup(chrom.tobytes(), stop_pos - 2, stop_pos, max_depth=100000000):

            if pileupcolumn.reference_pos == stop_pos - 1:
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
                            mut_tags.append(tag)
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
    bam.close()

    with open(json_file, "wb") as f:
        json.dump((tag_dict, cvrg_dict), f)

    # create fastq from aligned reads
    with open(outfile, 'w') as out:
        with open(file3, 'r') as families:
            for line in families:
                line = line.rstrip('\n')
                splits = line.split('\t')
                tag = splits[0]

                if tag in mut_tags:
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
