#!/usr/bin/env python

from __future__ import print_function

import argparse
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', type=argparse.FileType())
    parser.add_argument('output_file', type=argparse.FileType('w'), nargs='?', default=sys.stdout)
    args = parser.parse_args()
    lines = args.input_file.readlines()
    i = len(lines) - 1
    trimmed = False
    # step backwards through the lines, removing all As until we find a non-A nucleotide
    while not trimmed:
        line = lines[i].upper().rstrip()
        for j in range(len(line) - 1, -1, -1):
            # walk backwards through the line, checking for a non-A (and non-space) character
            if line[j] not in ['A', ' ']:
                lines[i] = line[:j + 1] + '\n'
                trimmed = True
                break
        else:
            # we processed the whole line - all As - so we don't include this line in the output
            i -= 1
    args.output_file.write(''.join(lines[:i + 1]))
