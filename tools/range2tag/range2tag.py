"""range2tag.py

Author -- Gundula Povysil
Contact -- povysil@bioinf.jku.at

Takes a SAM file, start and stop positions as input and prints all tags
of reads that overlap with regions to user specified output file.
=======  ==========  =================  ================================
Version  Date        Author             Description
0.0.2    2018-05-15  Gundula Povysil    -
=======  ==========  =================  ================================

USAGE: python range2tag.py inputFile.sam ranges.txt outputFile.txt
"""

import argparse
import os
import re
import sys

import numpy as np


def make_argparser():
    parser = argparse.ArgumentParser(description='Takes a SAM file, start and stop positions as input and prints all tags of reads that overlap with regions to user specified output file.')
    parser.add_argument('inputFile',
                        help='SAM file with aligned reads.')
    parser.add_argument('rangesFile', default=None,
                        help='BED file with cromosome, start and stop positions.')
    parser.add_argument('outputFile',
                        help='Output TXT file with tags that are within specified regions.')
    return parser


def range2tag(argv):
    parser = make_argparser()
    args = parser.parse_args(argv[1:])
    inputFile = args.inputFile
    rangesFile = args.rangesFile
    outputFile = args.outputFile

    if os.path.isfile(inputFile) is False:
        print("Error: Could not find '{}'".format(inputFile))
        exit(0)

    if rangesFile != str(None):
        with open(rangesFile, 'r') as regs:
            range_array = np.genfromtxt(regs, skip_header=0, delimiter='\t', comments='#', dtype='string')
        print(range_array.ndim)

        if range_array.ndim == 0:
            print("Error: file has 0 lines")
            exit(2)

        if range_array.ndim == 1:
            chrList = range_array[0]
            start_posList = range_array[1].astype(int)
            stop_posList = range_array[2].astype(int)
            chrList = [chrList.tolist()]
            start_posList = [start_posList.tolist()]
            stop_posList = [stop_posList.tolist()]
        else:
            chrList = range_array[:, 0]
            start_posList = range_array[:, 1].astype(int)
            stop_posList = range_array[:, 2].astype(int)

        if len(start_posList) != len(stop_posList):
            print("start_positions and end_positions do not have the same length")
            exit(3)

    with open(inputFile, 'r') as sam:
        data_array = np.genfromtxt(sam, skip_header=0, delimiter='\t', usecols=range(10), comments='#', dtype='string')

    tags = np.array(data_array[:, 0])
    print(len(tags))
    if re.search('_', tags[0]):
        tags = [re.split('_', x)[0] for x in tags]
    ref_pos = np.array(data_array[:, 3]).astype(int)
    cigar = np.array(data_array[:, 5]).astype(str)
    ref_chr = np.array(data_array[:, 2]).astype(str)
    # ref_chr_next = np.array(data_array[:, 6]).astype(str)

    lst = []
    ind = []
    # ref_name_next = []
    if rangesFile != str(None):
        chrList = np.array(chrList)
        start_posList = np.array(start_posList).astype(int)
        stop_posList = np.array(stop_posList).astype(int)

        for chr, start_pos, stop_pos in zip(chrList, start_posList, stop_posList):
            start_pos = start_pos - 3
            stop_pos = stop_pos + 3
            mut_tags = None
            for t in range(0, len(tags)):
                if cigar[t] != "*":
                    if ref_chr[t] == chr:
                        c_split = re.split('([A-Z])', cigar[
                            t])
                        cigar_long = None

                        for i in range(1, len(c_split),
                                       2):
                            if cigar_long is None:
                                cigar_long = np.repeat(c_split[i], c_split[
                                    i - 1])
                            else:
                                cigar_long = np.concatenate((cigar_long, np.repeat(c_split[i], c_split[i - 1])),
                                                            axis=0)

                        pos = ref_pos[t]
                        if pos < stop_pos:
                            for j in range(0, len(cigar_long)):
                                if pos >= stop_pos:
                                    break
                                if cigar_long[j] in ("M", "D", "N"):
                                    pos += 1
                            if pos > start_pos:
                                if mut_tags is None:
                                    mut_tags = np.array((tags[t]))
                                else:
                                    mut_tags = np.vstack((mut_tags, np.array(tags[t])))
                                # ref_name_next.append(ref_chr_next[t])
            if mut_tags is not None:
                index = np.repeat("{}_{}_{}".format(chr, start_pos + 3, stop_pos - 3), mut_tags.size)
                ind.append(index)
                if mut_tags.size == 1:
                    lst.append(np.vstack([mut_tags.tolist()]))
                else:
                    lst.append(mut_tags)
            else:
                print("No tags found in region {}_{}_{}!".format(chr, start_pos + 3, stop_pos - 3))
    else:
        tags = np.array(tags)[np.where((cigar != "*"))[0]]
        ref_chr = np.array(ref_chr)[np.where((cigar != "*"))[0]]
        # ref_chr_next = np.array(ref_chr_next)[np.where((cigar != "*"))[0]]
        # only_aligned_tags = np.where((ref_chr != "*") | (ref_chr_next != "*"))[0]
        only_aligned_tags = np.where(ref_chr != "*")[0]
        mut_tags = np.array(np.array(tags)[only_aligned_tags])
        index = np.array(np.array(ref_chr)[only_aligned_tags])
        # ref_name_next.append(np.array(ref_chr_next)[only_aligned_tags])
        ind.append(index)
        lst.append(mut_tags)
        # ref_name_next = np.concatenate((ref_name_next))
    index = np.concatenate((ind))
    tags = np.concatenate((lst))
    # mut_tags = np.column_stack((index, tags, ref_name_next))
    mut_tags = np.column_stack((index, tags))
    np.savetxt(outputFile, mut_tags, fmt="%s")
    print(len(mut_tags))
    print("File saved under {} in {}!".format(outputFile, os.getcwd()))


if __name__ == '__main__':
    sys.exit(range2tag(sys.argv))
