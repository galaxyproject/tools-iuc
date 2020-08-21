#!/usr/bin/env python
######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################
from __future__ import print_function
import sys
from argparse import ArgumentParser
import pandas as pd
from scipy.stats import gmean


def generate_MFI(input_file_name, output_file_name, mfi_calc):
    flock_df = pd.read_table(input_file_name)
    if mfi_calc == "mfi":
        MFIs = flock_df.groupby('Population').mean().round(decimals=2)
    elif mfi_calc == "gmfi":
        MFIs = flock_df.groupby('Population').agg(lambda x: gmean(list(x))).round(decimals=2)
    else:
        MFIs = flock_df.groupby('Population').median().round(decimals=2)

    with open(output_file_name, "w") as outf:
        MFIs.to_csv(outf, sep="\t", float_format='%.0f')
    return


if __name__ == "__main__":
    parser = ArgumentParser(
             prog="removeColumns",
             description="Generate MFI from Flow Result file.")

    parser.add_argument(
            '-i',
            dest="input_file",
            required=True,
            help="File location for the Flow Result file.")

    parser.add_argument(
            '-M',
            dest="mfi_calc",
            required=True,
            help="what to calculate for centroids.")

    parser.add_argument(
            '-o',
            dest="output_file",
            required=True,
            help="File location for the MFI output file.")

    args = parser.parse_args()
    generate_MFI(args.input_file, args.output_file, args.mfi_calc)
