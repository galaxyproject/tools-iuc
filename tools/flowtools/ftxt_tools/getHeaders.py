#!/usr/bin/env python
######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################
#
# Version 1.2 - May 2018
# added leeway for files with different nb of headers
#


from __future__ import print_function
import sys

from argparse import ArgumentParser


def print_headers(files, filenames, outfile):
    header_table = {}
    for i, eachfile in enumerate(files):
        with open(eachfile, "r") as ef:
            headers = ef.readline().strip()
            header_table[filenames[i]] = headers.split("\t")

    h = 0
    for f in header_table:
        j = len(header_table[f]) + 1
        if j > h:
            h = j

    idx = [str(x) for x in range(1, h)]

    with open(outfile, "w") as outf:
        outf.write("Index\t")
        outf.write("\t".join(idx) + "\n")
        for f in header_table:
            if len(header_table[f]) < h:
                for k in range(len(header_table[f]), h-1):
                    header_table[f].append("")
                sys.stderr.write(str(len(header_table[f])))
            outf.write(f + "\t")
            outf.write("\t".join(header_table[f]) + "\n")
    return


if __name__ == "__main__":
    parser = ArgumentParser(
             prog="GetHeaders",
             description="Gets the headers of all files in given set.")

    parser.add_argument(
            '-i',
            dest="input_files",
            required=True,
            action='append',
            help="File location for the text files.")

    parser.add_argument(
            '-n',
            dest="file_names",
            required=True,
            action='append',
            help="File names.")

    parser.add_argument(
            '-o',
            dest="output_file",
            required=True,
            help="Name of the output file.")

    args = parser.parse_args()
    input_files = [f for f in args.input_files]
    file_names = [fn for fn in args.file_names]
    print_headers(input_files, file_names, args.output_file)
