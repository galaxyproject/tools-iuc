#!/usr/bin/env python

######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################

# version 1.1 -- August 2017
# added checks for consistency between input files
# and upper limit on nb of cluster to look at

from __future__ import print_function
import sys
import os
import logging
import fileinput
import pandas as pd
from argparse import ArgumentParser
from jinja2 import Environment, FileSystemLoader
from shutil import copyfile
from collections import defaultdict


def check_pops(mfi_file, stat1):
    df = pd.read_table(mfi_file)
    df1 = pd.read_table(stat1)
    nb_pop = len(set(df.Population))
    nb_pop1 = len(df1.columns) - 2
    if (nb_pop > 40):
        sys.stderr.write("There are " + str(nb_pop) + " in the input file.")
        sys.exit(1)
    if (nb_pop != nb_pop1):
        sys.exit(2)


def panel_to_json_string(panel):
    # from http://stackoverflow.com/questions/28078118/merge-many-json-strings-with-python-pandas-inputs
    def __merge_stream(key, stream):
        return '"' + key + '"' + ': ' + stream + ', '
    try:
        stream = '{'
        for item in panel.items:
            stream += __merge_stream(item, panel.loc[item, :, :].to_json())
        # take out extra last comma
        stream = stream[:-2]
        # add the final paren
        stream += '}'
    except:
        logging.exception('Panel Encoding did not work')
    return stream


def get_outliers(group, upper, lower):
    cat = group.name
    out = {}
    for marker in group:
        out[marker] = group[(group[marker] > upper.loc[cat][marker]) | (group[marker] < lower.loc[cat][marker])][marker]
    return out


def get_boxplot_stats(all_data, mfi_file, output_json):
    # modified code from http://bokeh.pydata.org/en/latest/docs/gallery/boxplot.html
    # Get initial MFI values
    mfi = pd.read_table(mfi_file)
    mfi = mfi.set_index('Population')

    df = pd.read_table(all_data)
    # check if ever some pops not in cs_files
    missing_pop = [x for x in mfi.index if x not in set(df.Population)]

    if (missing_pop):
        zeros = {}
        for m in df.columns:
            zeros[m] = [0 for x in missing_pop]
        tmpdf = pd.DataFrame(zeros)
        tmpdf.Population = missing_pop
        df = df.append(tmpdf)

    pops = df.groupby('Population')
    q1 = pops.quantile(q=0.25)
    q2 = pops.quantile(q=0.5)
    q3 = pops.quantile(q=0.75)
    iqr = q3 - q1
    upper = q3 + 1.5*iqr
    lower = q1 - 1.5*iqr
    resampled = False
    # get outliers
    out = pops.apply(get_outliers, upper, lower).dropna()
    outliers = defaultdict(dict)
    for population in set(df.Population):
        for marker in df.columns:
            if marker != 'Population':
                tmp_outliers = list(out[population][marker])
                if (len(list(out[population][marker])) > 100):
                    tmp_outliers = list(out[population][marker].sample(n=100))
                    resampled = True
                outliers[population][marker] = tmp_outliers
    outdf = pd.DataFrame(outliers)

    data = {'q1': q1,
            'q2': q2,
            'q3': q3,
            'upper': upper,
            'lower': lower,
            'outliers': outdf.T,
            'mfi': mfi}
    wp = pd.Panel(data)

    with open(output_json, "w") as js_all:
        js_all.write(panel_to_json_string(wp))

    return resampled


def cs_overview(input_file, input_mfi, init_mfi, output_file, output_dir, tools_dir, cs_files):
    os.mkdir(output_dir)

    env = Environment(loader=FileSystemLoader(tools_dir + "/templates"))
    template = env.get_template("csOverview.template")

    real_directory = output_dir.replace("/job_working_directory", "")
    context = {'outputDirectory': real_directory}
    overview = template.render(**context)
    with open(output_file, "w") as outf:
        outf.write(overview)

    cs_overview_file = output_dir + "/csOverview.tsv"
    copyfile(input_file, cs_overview_file)

    cs_overview_mfis = output_dir + "/csAllMFIs.tsv"
    copyfile(input_mfi, cs_overview_mfis)

    # Get all the data to calculate quantiles, IRC and outliers.
    tmp_all_data = "csAllData.tsv"
    with open(tmp_all_data, "a") as alldata:
        # assumes that the files have ran through flock and CS and therefore have the same headers
        df1 = pd.read_table(cs_files[0])
        df1.to_csv(alldata, sep="\t", header=True, index=False)
        for i in range(1, len(cs_files)):
            df = pd.read_table(cs_files[i])
            df.to_csv(alldata, sep="\t", header=False, index=False)

    cs_boxplot_data = output_dir + "/csBoxplotData.json"
    resampled = get_boxplot_stats(tmp_all_data, init_mfi, cs_boxplot_data)
    if resampled:
        to_find = '<div id="outlierWarning" style="display:none;">'
        to_replace = '<div id="outlierWarning">'
        ## yay python 2.7
        ro = fileinput.input(output_file, inplace=True, backup=".bak")
        for roline in ro:
            print(roline.replace(to_find, to_replace), end='')
        ro.close()

    return


if __name__ == "__main__":
    parser = ArgumentParser(
             prog="csOverview",
             description="Generate an overview plot of crossSample results.")

    parser.add_argument(
            '-i',
            dest="input_file",
            required=True,
            help="File location for the summary statistics from CrossSample.")

    parser.add_argument(
            '-I',
            dest="input_mfi",
            required=True,
            help="File location for the MFI summary statistics from CrossSample.")

    parser.add_argument(
            '-s',
            dest="cs_outputs",
            required=True,
            action='append',
            help="File location for the CrossSample output files.")

    parser.add_argument(
            '-o',
            dest="output_file",
            required=True,
            help="File location for the HTML output file.")

    parser.add_argument(
            '-m',
            dest="mfi",
            required=True,
            help="File location for the MFI from FLOCK.")

    parser.add_argument(
            '-d',
            dest="output_directory",
            required=True,
            help="Directory location for the html supporting files.")

    parser.add_argument(
            '-t',
            dest="tool_directory",
            required=True,
            help="Location of the Tool Directory.")

    args = parser.parse_args()

    cs_files = [f for f in args.cs_outputs]
    check_pops(args.mfi, args.input_file)
    cs_overview(args.input_file, args.input_mfi, args.mfi, args.output_file, args.output_directory, args.tool_directory, cs_files)
