#!/usr/bin/env python

# adapted from https://github.com/bgruening/galaxytools/blob/f96142ca5acea989b828d6c92940172355b7f125/tools/jbrowse2/maf2bed.py
# which was "painfully converted from b0rken perl from:"
# https://unpkg.com/browse/jbrowse-plugin-mafviewer@1.0.6/dist/

import argparse
import sys


def maf2bed(assembly_name, input, output):
    id = 0
    buffer = ''
    start = 0
    end = 0
    score = 0
    chrom = ''

    db = "%s." % assembly_name
    # Read input from stdin
    for line in input:
        line = line.strip()
        if not line:
            continue

        line = line.split()
        if line[0] == 's' and line[1].startswith(db):
            chrom = line[1]
            chrom = chrom.replace(db, '')
            start = int(line[2])
            end = int(line[2]) + int(line[3])
            line = line[1:]
            line = ':'.join(line)
            temp = line
            buffer = temp if buffer == '' else f"{buffer},{temp}"
        elif line[0] == 'a':
            score = int(line[1].split('=')[1])
            if id > 0:
                output.write('\t'.join([chrom, '%d' % start, '%d' % end, f"{assembly_name}_{id}", '%d' % score, buffer]) + '\n')
            id += 1
            buffer = ''
        elif line[0] == 's':
            line = line[1:]
            line = ':'.join(line)
            temp = line
            buffer = temp if buffer == '' else f"{buffer},{temp}"

    output.write('\t'.join([chrom, '%d' % start, '%d' % end, f"{assembly_name}_{id}", '%d' % score, buffer]) + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="", epilog="")
    parser.add_argument("assembly_name", help="Assembly name")
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('output', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()

    maf2bed(args.assembly_name, args.input, args.output)
