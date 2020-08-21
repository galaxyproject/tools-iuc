#!/usr/bin/env python
######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################
# version 2
from __future__ import print_function

import sys
import os
from argparse import ArgumentParser
import pandas as pd
from scipy.stats import gmean


def run_FLOCK(input_file, method, bins, density, output_file, mfi_file,
              mfi_calc, profile, tool_directory):
    run_command = tool_directory + "/bin/" + method + " " + input_file
    if bins:
        run_command += " " + bins
    if density:
        run_command += " " + density

    os.system(run_command)

    move_command = "mv flock_results.txt " + output_file
    os.system(move_command)

    # Here add some way to calculate the count and tack it on to profile file.
    flockdf = pd.read_table(output_file)
    if mfi_calc == "mfi":
        MFIs = flockdf.groupby('Population').mean().round(decimals=2)
    elif mfi_calc == "gmfi":
        MFIs = flockdf.groupby('Population').agg(lambda x: gmean(list(x))).round(decimals=2)
    else:
        MFIs = flockdf.groupby('Population').median().round(decimals=2)

    with open(mfi_file, "w") as outf:
        MFIs.to_csv(outf, sep="\t", float_format='%.0f')

    (events, columns) = flockdf.shape
    fstats = {}
    fstats['population'] = flockdf.iloc[:, -1:].iloc[:, 0]
    fstats['population_freq'] = fstats['population'].value_counts()
    fstats['population_freq_sort'] = fstats['population_freq'].sort_index()
    fstats['population_per'] = (fstats['population'].value_counts(normalize=True) * 100).round(decimals=2)
    fstats['population_per_sort'] = fstats['population_per'].sort_index()
    fstats['population_all'] = pd.concat([fstats['population_freq_sort'], fstats['population_per_sort']], axis=1)
    fstats['population_all'].columns = ['Count', 'Percentage']
    fstats['population_all']['Population_ID'] = fstats['population_all'].index

    flock_profile = pd.read_table('profile.txt')
    profile_pop = flock_profile.merge(fstats['population_all'], on='Population_ID')
    profile_pop.to_csv(profile, sep="\t", float_format='%.2f', index=False)

#    get_profile = "mv profile.txt " + profile
#    os.system(get_profile)
    return


if __name__ == "__main__":
    parser = ArgumentParser(
             prog="runFlockMFI",
             description="Run Flock on text file and generate centroid file")

    parser.add_argument(
            '-i',
            dest="input_file",
            required=True,
            help="File location for the FCS file.")

    parser.add_argument(
            '-m',
            dest="method",
            required=True,
            help="Run flock1 or flock2.")

    parser.add_argument(
            '-M',
            dest="mfi_calc",
            required=True,
            help="what to calculate for centroids.")

    parser.add_argument(
            '-b',
            dest="bins",
            required=False,
            help="Number of Bins.")

    parser.add_argument(
            '-d',
            dest="density",
            required=False,
            help="Density.")

    parser.add_argument(
            '-o',
            dest="output_file",
            required=True,
            help="File location for the output file.")

    parser.add_argument(
            '-t',
            dest="tool_directory",
            required=True,
            help="File location for the output file.")

    parser.add_argument(
            '-c',
            dest="centroids",
            required=True,
            help="File location for the output centroid file.")

    parser.add_argument(
            '-p',
            dest="profile",
            required=True,
            help="File location for the output profile file.")

    args = parser.parse_args()
    run_FLOCK(args.input_file, args.method, args.bins,
              args.density, args.output_file, args.centroids, args.mfi_calc,
              args.profile, args.tool_directory)
