#!/usr/bin/env python

import argparse
import re


def padfile(infile, outfile, fieldcnt=None):
    with open(infile, 'r') as fh:
        out = open(outfile, 'w')
        commentlines = []
        tabs = '\t' * fieldcnt if fieldcnt is not None else None
        def pad_line(txtline, tabs=None):
            line = txtline.rstrip('\r\n')
            fields = line.split('\t')
            if not tabs:
                tabs = '\t' * len(fields)
            out.write('%s%s\n' % (line, tabs[len(fields):]))
        for i, txtline in enumerate(fh):
            if txtline.lstrip().startswith('#'):
                commentlines.append(txtline)
            else:
                if commentlines:
                    for i in range(len(commentlines)-1):
                        out.write(commentlines[i])
                    pad_line(commentlines[-1], tabs=tabs)
                    commentlines = []
                pad_line(txtline, tabs=tabs)
        out.close()


def fieldcount(infile):
    fieldcnt = 0
    with open(infile, 'r') as fh:
        for i, line in enumerate(fh):
            fieldcnt = max(fieldcnt, len(line.rstrip('\r\n').split('\t')))
    return fieldcnt


def tsvname(infile):
    return re.sub('\.txt$', '', infile) + '.tsv'


def __main__():
    parser = argparse.ArgumentParser(
        description='Pad a file with TABS for equal field size across lines')
    parser.add_argument(
        '-i', '--input', help='input file')
    parser.add_argument(
        '-o', '--output', help='output file')
    parser.add_argument(
        'files', nargs='*', help='.txt files')
    args = parser.parse_args()

    if args.input:
        outfile = args.output if args.output else tsvname(args.input)
        fieldcnt = fieldcount(args.input)
        padfile(args.input, outfile, fieldcnt=fieldcnt)
    for infile in args.files:
        outfile = tsvname(infile)
        fieldcnt = fieldcount(infile)
        padfile(infile, outfile, fieldcnt=fieldcnt)


if __name__ == "__main__":
    __main__()
