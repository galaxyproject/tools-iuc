#!/usr/bin/env python

"""
Input: Glimmer3 prediction
Output: GFF3 file
Return a GFF3 file with the genes predicted by Glimmer3
Bjoern Gruening

Note: Its not a full-fledged GFF3 file, its a really simple one.

"""

import re
import sys


def __main__():
    input_file = open(sys.argv[1], 'r')

    print('##gff-version 3\n')
    for line in input_file:
        line = line.strip()
        if line[0] == '>':
            header = line[1:]
        else:
            (id, start, end, frame, score) = re.split(r'\s+', line)
            if int(end) > int(start):
                strand = '+'
            else:
                strand = '-'
                (start, end) = (end, start)

            rest = 'frame=%s;score=%s' % (frame, score)
            print('\t'.join([header, 'glimmer_prediction', 'predicted_gene', start, end, '.', strand, '.', rest]))


if __name__ == "__main__":
    __main__()
