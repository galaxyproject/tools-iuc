#!/usr/bin/env python
######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################
from __future__ import print_function
import sys
import os
from scipy.stats import gmean
from argparse import ArgumentParser
from collections import defaultdict
import pandas as pd

#
# version 1.1 -- April 2016 -- C. Thomas
# modified to read in several input files and output to a directory
# + generates summary statistics
# also checks before running that input files are consistent with centroid file
#


def compare_MFIs(input_files, f_names, mfi_file):
    header_MFIs = ""
    flag_error = False
    with open(mfi_file, "r") as mfi_check:
        mfi_fl = mfi_check.readline().split("\t")
        header_MFIs = "\t".join([mfi_fl[h] for h in range(1, len(mfi_fl))])

    for hh, files in enumerate(input_files):
        with open(files, "r") as inf:
            hdrs = inf.readline()
            if hdrs != header_MFIs:
                sys.stderr.write(hdrs + "headers in " + f_names[hh] + " are not consistent with FLOCK centroid file:\n" + header_MFIs + "\n")
                flag_error = True
    if flag_error:
        sys.exit(2)


def stats_MFIs(cs_df, ctr, mfi_calc):
    if mfi_calc == "mfi":
        MFIs = cs_df.groupby('Population').mean().round(decimals=2)
    elif mfi_calc == "gmfi":
        MFIs = cs_df.groupby('Population').agg(lambda x: gmean(list(x))).round(decimals=2)
    else:
        MFIs = cs_df.groupby('Population').median().round(decimals=2)
    pop_freq = (cs_df.Population.value_counts(normalize=True) * 100).round(decimals=2)
    sorted_pop_freq = pop_freq.sort_index()
    MFIs['Percentage'] = sorted_pop_freq
    MFIs['Population'] = MFIs.index
    MFIs['SampleName'] = "".join(["Sample", str(ctr).zfill(2)])
    return MFIs


def get_pop_prop(input_files, summary_stat, mfi_stats, marker_names, mfi_calc):
    pop_count = defaultdict(dict)
    mrk = marker_names.strip().split("\t")
    markers = "\t".join([mrk[m] for m in range(1, len(mrk))])

    ctr_mfi = 0
    nb_pop = 0
    tot = {}
    with open(mfi_stats, "a") as mfis:
        mfis.write("\t".join([markers, "Percentage", "Population", "SampleName"]) + "\n")
        for files in input_files:
            cs = pd.read_table(files)
            tot[files] = len(cs.index)
            for pops in cs.Population:
                if pops in pop_count[files]:
                    pop_count[files][pops] += 1
                else:
                    pop_count[files][pops] = 1
            max_nb_pop = max(set(cs.Population))
            if (max_nb_pop > nb_pop):
                nb_pop = max_nb_pop
            ctr_mfi += 1
            cs_stats = stats_MFIs(cs, ctr_mfi, mfi_calc)
            cs_stats.to_csv(mfis, sep="\t", header=False, index=False)

    ctr = 0
    with open(summary_stat, "w") as outf:
        itpop = [str(x) for x in range(1, nb_pop + 1)]
        cols = "\t".join(itpop)
        outf.write("FileID\tSampleName\t" + cols + "\n")
        for eachfile in pop_count:
            tmp = []
            for num in range(1, nb_pop + 1):
                if num not in pop_count[eachfile]:
                    pop_count[eachfile][num] = 0
                tmp.append(str((pop_count[eachfile][num] / float(tot[eachfile])) * 100))
            props = "\t".join(tmp)
            ctr += 1
            sample_name = "".join(["Sample", str(ctr).zfill(2)])
            outf.write("\t".join([input_files[eachfile], sample_name, props]) + "\n")


def run_cross_sample(input_files, f_names, mfi_file, output_dir, summary_stat,
                     mfi_stats, tool_directory, mfi_calc):
    markers = ""
    # Strip off Header Line
    with open(mfi_file, "r") as mfi_in, open("mfi.txt", "w") as mfi_out:
        markers = mfi_in.readline().strip("\n")
        for line in mfi_in:
            mfi_out.write(line)

    # Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    outputs = {}
    # Run cent_adjust
    for nm, flow_file in enumerate(input_files):
        run_command = tool_directory + "/bin/cent_adjust mfi.txt " + flow_file
        print(run_command)
        os.system(run_command)
        flow_name = os.path.split(flow_file)[1]
        outfile = os.path.join(output_dir, flow_name + ".flowclr")
        outputs[outfile] = f_names[nm]
        with open(flow_file, "r") as flowf, open("population_id.txt", "r") as popf, open(outfile, "w") as outf:
            f_line = flowf.readline()
            f_line = f_line.rstrip()
            f_line = f_line + "\tPopulation\n"
            outf.write(f_line)

            for line in flowf:
                line = line.rstrip()
                pop_line = popf.readline()
                pop_line = pop_line.rstrip()
                line = line + "\t" + pop_line + "\n"
                outf.write(line)
    get_pop_prop(outputs, summary_stat, mfi_stats, markers, mfi_calc)
    return


def generate_CS_stats(mfi_stats, all_stats):
    df = pd.read_table(mfi_stats)
    means = df.groupby('Population').mean().round(decimals=2)
    medians = df.groupby('Population').median().round(decimals=2)
    stdev = df.groupby('Population').std().round(decimals=2)
    all_markers = []
    with open(mfi_stats, "r") as ms:
        ms_fl = ms.readline().strip()
        all_markers = ms_fl.split("\t")[0:-2]

    with open(all_stats, "w") as mstats:
        hdgs = ["\t".join(["_".join([mrs, "mean"]), "_".join([mrs, "median"]), "_".join([mrs, "stdev"])]) for mrs in all_markers]
        mstats.write("Population\t")
        mstats.write("\t".join(hdgs) + "\n")
        for pops in set(df.Population):
            tmp_line = []
            for mar in all_markers:
                tmp_line.append("\t".join([str(means.loc[pops, mar]), str(medians.loc[pops, mar]), str(stdev.loc[pops, mar])]))
            mstats.write(str(pops) + "\t")
            mstats.write("\t".join(tmp_line) + "\n")


if __name__ == "__main__":
    parser = ArgumentParser(
             prog="runCrossSample",
             description="Run CrossSample on Flow file")

    parser.add_argument(
            '-i',
            dest="input_files",
            required=True,
            action='append',
            help="File locations for flow text files.")

    parser.add_argument(
            '-n',
            dest="filenames",
            required=True,
            action='append',
            help="Filenames")

    parser.add_argument(
            '-m',
            dest="mfi",
            required=True,
            help="File location for the MFI text file.")

    parser.add_argument(
            '-o',
            dest="out_path",
            required=True,
            help="Path to the directory for the output files.")

    parser.add_argument(
            '-M',
            dest="mfi_calc",
            required=True,
            help="what to calculate for centroids.")

    parser.add_argument(
            '-s',
            dest="sstat",
            required=True,
            help="File location for the summary statistics.")

    parser.add_argument(
            '-S',
            dest="mfi_stat",
            required=True,
            help="File location for the MFI summary statistics.")

    parser.add_argument(
            '-t',
            dest="tool_dir",
            required=True,
            help="File location for cent_adjust.")

    parser.add_argument(
            '-a',
            dest="all_stats",
            required=True,
            help="File location for stats on all markers.")

    args = parser.parse_args()

    input_files = [f for f in args.input_files]
    input_names = [n for n in args.filenames]
    compare_MFIs(input_files, input_names, args.mfi)
    run_cross_sample(input_files, input_names, args.mfi, args.out_path, args.sstat, args.mfi_stat, args.tool_dir, args.mfi_calc)
    generate_CS_stats(args.mfi_stat, args.all_stats)
