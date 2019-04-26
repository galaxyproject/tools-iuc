#!/usr/bin/env python

######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################

from __future__ import print_function
from __future__ import division
import sys
import os
import pandas as pd
from argparse import ArgumentParser


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def compare_headers(files):
    headers = {}
    for eachfile in files:
        with open(eachfile, "r") as ef:
            headers[eachfile] = ef.readline().strip().lower().split("\t")

    hdgs_in_common = []
    flag = {}

    for ref_hdgs in headers[files[0]]:
        flag[ref_hdgs] = 1

        for ij in range(1, len(files)):
            if ref_hdgs in headers[files[ij]]:
                flag[ref_hdgs] += 1
        if flag[ref_hdgs] == len(files):
            hdgs_in_common.append(ref_hdgs)

    if not hdgs_in_common:
        sys.exit(9)
    return(hdgs_in_common)


def get_nb_lines(files):
    tot_event = 0
    for f in files:
        df = pd.read_table(f)
        tot_event += (len(df.index) - 1)
    return(tot_event)


def get_headers_index(list_headings, headings):
    idxs = []
    lhdgs = [x.lower() for x in headings]
    for element in list_headings:
        idxs.append(int(lhdgs.index(element)))
    return(idxs)


def merge_and_DS_txt(in_files, out_file, col_names, factor_ds):
    """Concatenates together tab-separated files.
    The output will have only the columns in common to all the files provided
    as input, as determined by the headers.
    All lines after the header line must contain only numbers.
    Potential errors are logged to stderr. If the number of errors reaches 10,
    the program stops.
    If a downsampling factor is given, returns the indicated fraction of
    random lines.
    """

    nb_errors = 0
    max_error = 10

    # get list of headers in common to all files
    list_hdgs = compare_headers(in_files)
    total_events = get_nb_lines(in_files)
    total_final = total_events * ds_factor
    nb_per_file = int(total_final / len(in_files))

    with open(out_file, "w") as outf:
        ff_order = []
        # HEADERS:
        with open(in_files[0], "r") as first_file:
            headings_ff = first_file.readline().strip()
            headings = headings_ff.split("\t")
            # Get index of headers in common:
            hdrs_idx = get_headers_index(list_hdgs, headings)

            # If column to merge on were provided:
            if col_names:
                for ix in col_names:
                    if ix not in hdrs_idx:
                        nb_errors += 1
                        sys.stderr.write(" ".join(["WARNING: column", str(ix), "in", in_files[0],
                                                   "does not exist in all files or has a different header.\n"]))
                        if nb_errors == max_error:
                            exit_code = 4
                            sys.stderr.write("Run aborted - too many errors.")
                            os.remove(out_file)
                hdrs_idx = col_names

            # Print out to output file:
            headings_to_write = []
            for cti in range(0, len(headings)):
                if cti in hdrs_idx:
                    headings_to_write.append(headings[cti])
                    ff_order.append(headings[cti])
            outf.write("\t".join(headings_to_write) + "\n")

        # DATA
        for infile in in_files:
            with open(infile, "r") as inf:
                headings_inf = inf.readline().strip()
                hdgs = headings_inf.split("\t")
                # Get the index of columns to keep:
                hdgs_idx = []
                for ctc in ff_order:
                    hdgs_idx.append(int(hdgs.index(ctc)))
                if col_names:
                    for iy in col_names:
                        if iy not in hdgs_idx:
                            nb_errors += 1
                            sys.stderr.write(" ".join(["WARNING: column", str(iy), "in", infile,
                                                       "does not exist in all files or has a different header.\n"]))
                            if nb_errors == max_error:
                                exit_code = 4
                                sys.stderr.write("Run aborted - too many errors.")
                                os.remove(out_file)
                    hdgs_idx = col_names

            df = pd.read_table(infile, usecols=hdrs_idx)
            df_ds = df.sample(nb_per_file, replace=False)

            for cols in df_ds.columns.values:
                if df_ds[cols].count() != len(df_ds[cols]):
                    sys.stderr.write(infile + "contains non-numeric data\n")

                    with open(infile, "r") as checkfile:
                        fl = checkfile.readline()
                        count_lines = 1
                        for checklines in checkfile:
                            to_check = checklines.strip().split("\t")
                            count_lines += 1
                            for item in to_check:
                                if not is_number(item):
                                    sys.stderr.write(" ".join(["WARNING: line", str(count_lines),
                                                               "in", infile, "contains non-numeric results\n"]))
                    sys.exit(2)

            df_ds = df_ds.ix[:, ff_order]
            df_ds.to_csv(outf, sep="\t", header=False, index=False)

    if nb_errors > 0:
        exit_code = 3
        if nb_errors == max_error:
            exit_code = 4
            sys.stderr.write("Run aborted - too many errors.")
            os.remove(out_file)
        sys.exit(exit_code)
    return


if __name__ == "__main__":
    parser = ArgumentParser(
             prog="FCStxtmerge",
             description="Merge based on headers text-converted FCS files into one text file.")

    parser.add_argument(
            '-i',
            dest="input_files",
            required=True,
            action='append',
            help="File location for the text files.")

    parser.add_argument(
            '-o',
            dest="output_file",
            required=True,
            help="Name of the output file.")

    parser.add_argument(
            '-c',
            dest="columns",
            help="Specify which column to keep in output file")

    parser.add_argument(
            '-d',
            dest="downsampling_factor",
            help="How much of each file to keep")

    args = parser.parse_args()

    # Get columns to merge on if any:
    default_value_col = ["i.e.:1,2,5", "default", "Default"]
    columns = []
    if args.columns:
        if args.columns not in default_value_col:
            tmp_col = args.columns.split(",")
            if len(tmp_col) == 1:
                if not tmp_col[0].strip():
                    columns = []
                elif not is_integer(tmp_col[0].strip()):
                    sys.exit(7)
                else:
                    columns.append(int(tmp_col[0].strip()) - 1)
            else:
                for c in range(0, len(tmp_col)):
                    if not is_integer(tmp_col[c].strip()):
                        sys.exit(6)
                    else:
                        columns.append(int(tmp_col[c].strip()) - 1)

    # Get down sampling factor if any:
    # Note: change '%' to 'X' because somehow that's what Galaxy passes?
    default_value_ds = ["i.e.:0.1 or 10X", "default", "Default"]
    ds_factor = 0.1
    if args.downsampling_factor:
        if args.downsampling_factor not in default_value_ds:
            args.downsampling_factor = args.downsampling_factor.strip()
            downsampling_factor = args.downsampling_factor.rstrip("X")
            if is_number(downsampling_factor):
                ds_factor = float(downsampling_factor)
                if ds_factor > 1 and ds_factor <= 100:
                    ds_factor = float(downsampling_factor) / 100
                elif ds_factor > 100 or ds_factor <= 0:
                    sys.stderr.write(str(ds_factor))
                    sys.exit(8)
            else:
                sys.exit(8)

    input_files = [f for f in args.input_files]
    merge_and_DS_txt(input_files, args.output_file, columns, ds_factor)
