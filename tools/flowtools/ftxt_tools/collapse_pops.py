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


def collapse_populations(in_file, out_file, populations, collapse_in, exit_code):
    df = pd.read_table(in_file, dtype={'Population': object})
    df['new_population'] = df.Population

    for i, sets_pop in enumerate(populations):
        df.loc[df['Population'].isin(sets_pop), ['new_population']] = collapse_in[i]

    df.Population = df.new_population
    df.drop(['new_population'], inplace=True, axis=1)

    df.to_csv(out_file, sep="\t", index=False)

    sys.exit(exit_code)


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
            action='append',
            help="List of populations to collapse.")

    parser.add_argument(
            '-o',
            dest="output_file",
            required=True,
            help="Name of the output file.")

    parser.add_argument(
            '-c',
            dest="collapse_pop",
            required=True,
            action='append',
            help="What to collapse the populations in.")

    args = parser.parse_args()

    # check populations
    default_values_pop = ["i.e.:2,3,11,25", "default", "Default"]
    default_values_col = ["i.e.:4", "default", "Default"]
    pops = [p for p in args.pops]
    popc = [pc.strip() for pc in args.collapse_pop]
    exit_code = 0
    # Check sets of pops to collapse
    populations = []
    total_pops = []
    for pop_set in pops:
        if pop_set not in default_values_pop:
            tmp_pops = pop_set.split(",")
            for popn in tmp_pops:
                if not is_int(popn):
                    sys.exit(6)
                else:
                    total_pops.append(int(popn))
            strp_pops = [p.strip() for p in tmp_pops]
            populations.append(strp_pops)
        else:
            sys.exit(4)
    if len(total_pops) != len(set(total_pops)):
        sys.exit(7)
    # Check pops to collapse in
    collapse_in = []
    for col_pop in popc:
        if col_pop not in default_values_col:
            if not is_int(col_pop):
                sys.exit(6)
            else:
                if int(col_pop) > 40:
                    exit_code = 2
                collapse_in.append(col_pop)
        else:
            sys.exit(6)
    if len(collapse_in) != len(set(collapse_in)):
        exit_code += 3

    collapse_populations(args.input_file, args.output_file, populations, collapse_in, exit_code)
