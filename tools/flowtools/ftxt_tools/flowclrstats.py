#!/usr/bin/env python

######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################

from __future__ import print_function
import sys
from argparse import ArgumentParser
import pandas as pd


def get_FLOCK_stats(input_file, output_file, out_file2):
    df = pd.read_table(input_file)
    summary = df.groupby('Population').describe().round(1)
    counts = df['Population'].value_counts()
    percent = (df['Population'].value_counts(normalize=True) * 100).round(decimals=2)
    tot_count = len(df['Population'])

    to_rm = summary.loc(axis=0)[:, ['count']].index.tolist()
    df1 = summary[~summary.index.isin(to_rm)]
    df1.to_csv(out_file2, sep="\t")

    with open(output_file, "w") as outf:
        outf.write("Population\tCount\tPercentage\n")
        for pops in set(df.Population):
            outf.write("\t".join([str(pops), str(counts.loc[pops]), str(percent.loc[pops])]) + "\n")
        outf.write("Total\t" + str(tot_count) + "\t \n")
    return


if __name__ == '__main__':
    parser = ArgumentParser(
            prog="flowstats",
            description="Gets statistics on FLOCK run")

    parser.add_argument(
          '-i',
          dest="input_file",
          required=True,
          help="File locations for flow clr file.")

    parser.add_argument(
          '-o',
          dest="out_file",
          required=True,
          help="Path to the directory for the output file.")

    parser.add_argument(
          '-p',
          dest="out_file2",
          required=True,
          help="Path to the directory for the output file.")
    args = parser.parse_args()

    get_FLOCK_stats(args.input_file, args.out_file, args.out_file2)
