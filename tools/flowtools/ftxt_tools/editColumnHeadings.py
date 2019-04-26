#!/usr/bin/env python

######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################
#
# Cristel Thomas - May 2018
# Version 2 -- with Pandas!
#

from __future__ import print_function
import sys

from argparse import ArgumentParser
import pandas as pd


def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def rearrange_file(input_file, output_file, new_cols, new_order, flag_text):
    df = pd.read_table(input_file)
    original_columns = [x for x in df.columns]
    if new_cols:
        edited_cols = []
        if len(new_cols) > len(df.columns):
            sys.exit(6)
        for i in range(0, len(df.columns)):
            if df.columns[i] in new_cols:
                edited_cols.append(new_cols[df.columns[i]])
            else:
                edited_cols.append(df.columns[i])
        df.columns = edited_cols

    if new_order:
        if len(new_order) > len(df.columns):
            sys.exit(6)
        subset = []
        if flag_text:
            check_subset = sum([False if x in df.columns else True for x in new_order])
            if check_subset > 0:
                sys.exit(9)
            subset = new_order
        else:
            subset = [df.columns[x] for x in new_order]
        df = df[subset]

    df.to_csv(output_file, sep="\t", index=False)
    if new_cols:
        for c in new_cols:
            if c not in original_columns:
                sys.exit(10)

if __name__ == "__main__":
    parser = ArgumentParser(
             prog="editColumnHeadings",
             description="Cut, rearrange and rename columns in a tab-separated file.")

    parser.add_argument(
            '-i',
            dest="input_file",
            required=True,
            help="File location for the text file.")

    parser.add_argument(
            '-r',
            dest="columns",
            action="append",
            help="Columns to replace.")

    parser.add_argument(
            '-w',
            dest="replace_with",
            action="append",
            help="new column headers.")

    parser.add_argument(
            '-n',
            dest="new_order",
            help="New column order if re-ordering or subsetting.")

    parser.add_argument(
            '-o',
            dest="output_file",
            required=True,
            help="Name of the output file.")

    args = parser.parse_args()


    new_order = []
    new_cols = {}
#    flag = False
#    exit_codes = [3,4,7,8,9,10,2]
    defaults = ["i.e.:TLR 6, TLR6PE", "i.e.:TLR6", "i.e.:1,2,5 or CD3,CD4,CCR3", "default", "Default", ""]
    flag_text = False

    if args.new_order:
        if args.new_order not in defaults:
            nwor = [x.strip() for x in args.new_order.strip().split(",")]
            check_integer = [is_integer(x) for x in nwor]
            if sum(check_integer) != len(check_integer):
                flag_text = True
            new_order = [str(x) if flag_text else int(x)-1 for x in nwor]
        else:
            sys.exit(8)

    if args.columns:
        if args.replace_with:
            cols_to_change = [c.strip().split(",") if c not in defaults else None for c in args.columns]
            replacements = [r.strip() if r not in defaults else None for r in args.replace_with]
            check_col = sum([True if x is not None else False for x in cols_to_change])
            check_rep = sum([True if x is not None else False for x in replacements])
            if check_col != check_rep:
                sys.exit(7)
            for i in range(0, check_col):
                if cols_to_change[i]:
                    if replacements[i]:
                        for c in cols_to_change[i]:
                            new_cols[c.strip()] = replacements[i]
                    else:
                        sys.exit(4)
                else:
                    sys.exit(3)
        else:
            sys.exit(7)
    else:
        if args.replace_with:
            sys.exit(7)

    if not new_order and not new_cols:
        sys.exit(2)

    rearrange_file(args.input_file, args.output_file, new_cols, new_order, flag_text)
