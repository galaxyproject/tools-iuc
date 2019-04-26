#!/usr/bin/env python
######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################
from __future__ import print_function
from __future__ import division
import pandas as pd
from argparse import ArgumentParser
import sys


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def error_report(input_file, fname, output_file):
    errors = 0
    df = pd.read_table(input_file)
    with open(output_file, "w") as outf:
        for cols in df.columns.values:
            if df[cols].count() != len(df[cols]):
                with open(input_file, "r") as checkfile:
                    fl = checkfile.readline()
                    count_lines = 1
                    for checklines in checkfile:
                        to_check = checklines.strip().split("\t")
                        count_lines += 1
                        for item in to_check:
                            if not is_number(item):
                                errors += 1
                                outf.write(" ".join(["WARNING: line", str(count_lines), "in", fname, "contains non-numeric results\n"]))
        if errors == 0:
            outf.write("No errors in the file.\n")
    return


if __name__ == "__main__":
    parser = ArgumentParser(
             prog="txtDiagnosis",
             description="Reports potential errors in text-converted FCS files")

    parser.add_argument(
            '-i',
            dest="input_file",
            required=True,
            help="File location for the text file.")

    parser.add_argument(
            '-n',
            dest="filename",
            required=True,
            help="Filename location for the text file.")

    parser.add_argument(
            '-o',
            dest="output_file",
            required=True,
            help="Name of the output file.")

    args = parser.parse_args()

    error_report(args.input_file, args.filename, args.output_file)
