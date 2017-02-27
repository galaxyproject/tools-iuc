#!/usr/bin/env python

######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################

from __future__ import print_function
import sys

from argparse import ArgumentParser


def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def rearrange_file(input_file, col_order, col_names, output_file):
    with open(input_file, "r") as infl, open(output_file, "w") as outf:
        # headers
        hdrs = infl.readline().strip()
        current_hdrs = hdrs.split("\t")
        if not col_order and col_names:
            if len(col_names) != len(current_hdrs):
                sys.stderr.write("There are " + str(len(current_hdrs)) + " columns but " + str(len(col_names)) + " marker names were provided\n")
                sys.exit(4)
        if col_names:
            tmp_hdr = []
            for i in range(0, len(col_names)):
                if col_names[i].strip():
                    tmp_hdr.append(col_names[i].strip())
                else:
                    if col_order:
                        tmp_hdr.append(current_hdrs[col_order[i]])
                    else:
                        tmp_hdr.append(current_hdrs[i])
            hdrs = ("\t".join(tmp_hdr))
        elif col_order:
            tp_hdr = []
            for j in col_order:
                tp_hdr.append(current_hdrs[j])
            hdrs = ("\t".join(tp_hdr))

        outf.write(hdrs + "\n")

        # columns
        for lines in infl:
            cols = lines.strip().split("\t")
            if not col_order:
                col_order = [x for x in range(0, len(current_hdrs))]
            outf.write("\t".join([cols[c] for c in col_order]) + "\n")


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
            '-c',
            dest="columns",
            help="Columns to keep in the order to keep them in.")

    parser.add_argument(
            '-n',
            dest="column_names",
            help="Column names if renaming.")

    parser.add_argument(
            '-o',
            dest="output_file",
            required=True,
            help="Name of the output file.")

    args = parser.parse_args()

    # check column indices
    default_value_col = ["i.e.:1,5,2", "default", "Default"]
    col_order = []
    if args.columns:
        if args.columns not in default_value_col:
            tmp_col = args.columns.split(",")
            if len(tmp_col) == 1:
                if not tmp_col[0].strip():
                    col_order = []
                elif not is_integer(tmp_col[0].strip()):
                    sys.exit(2)
                else:
                    col_order.append(int(tmp_col[0].strip()) - 1)
            else:
                for c in range(0, len(tmp_col)):
                    if not is_integer(tmp_col[c].strip()):
                        sys.exit(3)
                    else:
                        col_order.append(int(tmp_col[c].strip()) - 1)

    # check column names
    default_value_nms = ["i.e.:Marker1,,Marker4", "default", "Default"]
    col_names = []
    if args.column_names:
        if args.column_names not in default_value_nms:
            col_names = args.column_names.split(",")
            if col_order:
                if len(col_order) != len(col_names):
                    sys.stderr.write("There are " + str(len(col_order)) + " columns selected and " + str(len(col_names)) + " marker names\n")
                    sys.exit(4)

    rearrange_file(args.input_file, col_order, col_names, args.output_file)

    sys.exit(0)
