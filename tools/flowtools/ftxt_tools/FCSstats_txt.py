#!/usr/bin/env python

######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################

from __future__ import print_function
import pandas as pd
from argparse import ArgumentParser


def get_txt_stats(in_file, out_file):
    df = pd.read_table(in_file)
    summary = df.describe().round(1)
    df1 = summary[1:]
    x = summary[:1].values.tolist()
    df1.to_csv(out_file, sep="\t")
    with open(out_file, "a") as ot:
        ot.write("\n\n" + str(int(x[0][0])) + " events\n")
    return


if __name__ == "__main__":
    parser = ArgumentParser(
             prog="getTxtStats",
             description="Prints summary statistics from given file.")

    parser.add_argument(
            '-i',
            dest="input_file",
            required=True,
            help="File location for the text file.")

    parser.add_argument(
            '-o',
            dest="output_file",
            required=True,
            help="Name of the output file.")

    args = parser.parse_args()

    get_txt_stats(args.input_file, args.output_file)
