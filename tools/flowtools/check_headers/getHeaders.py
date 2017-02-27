#!/usr/bin/env python
######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################
from __future__ import print_function
import sys

from argparse import ArgumentParser


def print_headers(files, filenames, outfile):
    with open(outfile, "w") as outf:
        for i, eachfile in enumerate(files):
            with open(eachfile, "r") as ef:
                headers = ef.readline()
                outf.write("\t".join([filenames[i], headers]))
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
    sys.exit(0)
