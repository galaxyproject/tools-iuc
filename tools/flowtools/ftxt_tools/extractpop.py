#!/usr/bin/env python

######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################

from __future__ import print_function
import sys
import pandas as pd

from argparse import ArgumentParser


def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def extract_pop(in_file, pop_list, out_file):
    df = pd.read_table(in_file, dtype={'Population': object})
    dfout = df.loc[df['Population'].isin(pop_list)]
    dfout.to_csv(out_file, sep="\t", index=False)
    return


def remove_pop(in_file, pop_list, out_file):
    df = pd.read_table(in_file, dtype={'Population': object})
    dfout = df.loc[~df['Population'].isin(pop_list)]
    dfout.to_csv(out_file, sep="\t", index=False)
    return


if __name__ == "__main__":
    parser = ArgumentParser(
             prog="ExtractPop",
             description="Extract events associated to given population numbers.")

    parser.add_argument(
            '-i',
            dest="input_file",
            required=True,
            help="File location for the text file.")

    parser.add_argument(
            '-p',
            dest="pops",
            required=True,
            help="List of populations to extract.")

    parser.add_argument(
            '-o',
            dest="output_file",
            required=True,
            help="Name of the output file.")

    parser.add_argument(
            '-m',
            dest="method",
            required=True,
            help="What to do with the populations.")

    args = parser.parse_args()

    # check populations
    default_values = ["i.e.:2,3,11,25", "default", "Default"]
    populations = []
    if args.pops:
        if args.pops not in default_values:
            tmp_pops = args.pops.split(",")
            for popn in tmp_pops:
                populations.append(popn.strip())
        else:
            sys.exit(2)
    for pops in populations:
        if not is_int(pops):
            sys.exit(3)
    if args.method == "selected":
        extract_pop(args.input_file, populations, args.output_file)
    else:
        remove_pop(args.input_file, populations, args.output_file)
