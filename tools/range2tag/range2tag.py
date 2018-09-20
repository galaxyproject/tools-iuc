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
import re
import sys
import os

import numpy as np


def make_argparser():
    parser = argparse.ArgumentParser(description='Takes a SAM file, start and stop positions as input and prints all tags of reads that overlap with regions to user specified output file.')
    parser.add_argument('inputFile',
                        help='SAM file with aligned reads.')
    parser.add_argument('rangesFile',
                        help='TXT file with start and stop positions.')
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

    if os.path.isfile(rangesFile) is False:
        print("Error: Could not find '{}'".format(rangesFile))
        exit(0)

    with open(rangesFile, 'r') as regs:
        range_array = np.genfromtxt(regs, skip_header=0, delimiter='\t', comments='#')

    start_posList = range_array[:, 0].astype(int)
    stop_posList = range_array[:, 1].astype(int)

    if len(start_posList) == 0:
        print("Error: start_positions is empty")
        exit(2)

    if len(stop_posList) == 0:
        print("Error: end_positions is empty")
        exit(3)

    if len(start_posList) != len(stop_posList):
        print("start_positions and end_positions do not have the same length")
        exit(3)

    with open(inputFile, 'r') as sam:
        data_array = np.genfromtxt(sam, skip_header=0, delimiter='\t', usecols=range(11), comments='#', dtype='string')

    tags = np.array(data_array[:, 0])
    ref_pos = np.array(data_array[:, 3]).astype(int)
    cigar = np.array(data_array[:, 5])

    lst = []
    ind = []
    start_posList = np.array(start_posList).astype(int)
    stop_posList = np.array(stop_posList).astype(int)

    for start_pos, stop_pos in zip(start_posList, stop_posList):
        start_pos = start_pos - 3
        stop_pos = stop_pos + 3
        mut_tags = None
        for t in range(0, len(tags)):
            if cigar[t] != "*":
                c_split = re.split('([A-Z])', cigar[t])
                cigar_long = None

                for i in range(1, len(c_split), 2):
                    if cigar_long is None:
                        cigar_long = np.repeat(c_split[i], c_split[i - 1])
                    else:
                        cigar_long = np.concatenate((cigar_long, np.repeat(c_split[i], c_split[i - 1])), axis=0)

                pos = ref_pos[t]
                # seq_pos = 0
                #    print(pos)
                if pos < stop_pos:
                    for j in range(0, len(cigar_long)):
                        if pos >= stop_pos:
                            break
                        if cigar_long[j] in ("M", "D", "N"):
                            pos += 1
                            #        print(pos)
                    if pos > start_pos:
                        if mut_tags is None:
                            mut_tags = np.array((tags[t]))
                        else:
                            mut_tags = np.vstack((mut_tags, np.array(tags[t])))

        index = np.repeat("{}_{}".format(start_pos, stop_pos), len(mut_tags))
        ind.append(index)
        lst.append(mut_tags)

    index = np.concatenate((ind))
    tags = np.concatenate((lst))
    mut_tags = np.column_stack((index, tags))

    np.savetxt(outputFile, mut_tags, fmt="%s")
    print("File saved under {} in {}!".format(outputFile, os.getcwd()))


if __name__ == '__main__':
    sys.exit(range2tag(sys.argv))
