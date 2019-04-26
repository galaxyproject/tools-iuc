######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################
from __future__ import print_function
import sys
import pandas as pd
from scipy.stats import gmean
from argparse import ArgumentParser


def gen_overview_stats(file_name):
    flow_stats = {}
    fcs = pd.read_table(file_name)
    (events, columns) = fcs.shape
    flow_stats['fcs'] = fcs
    flow_stats['events'] = events
    flow_stats['columns'] = columns - 1
    flow_stats['data'] = fcs.iloc[:, :-1]
    flow_stats['population'] = fcs.iloc[:, -1:].iloc[:, 0]
    flow_stats['population_freq'] = flow_stats['population'].value_counts()
    flow_stats['population_sample'] = (flow_stats['population_freq'] * (20000/float(events))).round(decimals=0)
    flow_stats['population_freq_sort'] = flow_stats['population_freq'].sort_index()
    flow_stats['population_per'] = (flow_stats['population'].value_counts(normalize=True) * 100).round(decimals=2)
    flow_stats['population_per_sort'] = flow_stats['population_per'].sort_index()
    flow_stats['population_all'] = pd.concat([flow_stats['population_freq_sort'], flow_stats['population_per_sort']], axis=1)
    flow_stats['population_all'].columns = ['Count', 'Percentage']
    flow_stats['min'] = flow_stats['data'].values.min()
    flow_stats['max'] = flow_stats['data'].values.max()
    flow_stats['markers'] = list(flow_stats['data'].columns)
    flow_stats['mfi'] = fcs.groupby('Population').mean().round(decimals=2)
    flow_stats['mfi_pop'] = pd.merge(flow_stats['mfi'], flow_stats['population_all'], left_index=True, right_index=True)
    flow_stats['mfi_pop']['Population'] = flow_stats['mfi_pop'].index
    flow_stats['gmfi'] = fcs.groupby('Population').agg(lambda x: gmean(list(x))).round(decimals=2)
    flow_stats['gmfi_pop'] = pd.merge(flow_stats['gmfi'], flow_stats['population_all'], left_index=True, right_index=True)
    flow_stats['gmfi_pop']['Population'] = flow_stats['gmfi_pop'].index
    flow_stats['mdfi'] = fcs.groupby('Population').median().round(decimals=2)
    flow_stats['mdfi_pop'] = pd.merge(flow_stats['mdfi'], flow_stats['population_all'], left_index=True, right_index=True)
    flow_stats['mdfi_pop']['Population'] = flow_stats['mdfi_pop'].index

    #
    # If the number of events is less than 20000, then return
    # the complete data set,
    # Otherwise sample the data to only return 20000 events.
    if events <= 20000:
        flow_stats['sample'] = fcs
    else:
        fcs_np = fcs.values
        sample_data = []
        pop_found = {}
        for i in range(0, events):
            population_number = fcs_np[i][columns-1]
            if population_number in pop_found:
                if pop_found[population_number] < flow_stats['population_sample'][population_number]:
                    pop_found[population_number] += 1
                    sample_data.append(fcs_np[i])
            else:
                pop_found[population_number] = 1
                sample_data.append(fcs_np[i])
        flow_stats['sample'] = pd.DataFrame(sample_data)
        flow_stats['sample'].columns = fcs.columns

    flow_stats['sample_data'] = flow_stats['sample'].iloc[:, :-1]
    flow_stats['sample_population'] = flow_stats['sample'].iloc[:, -1:].iloc[:, 0]

    return flow_stats


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
    args = parser.parse_args()

    flow_stats = gen_overview_stats(args.input_file)
    with open(args.out_file, "w") as outf:
        outf.write("Events: ", flow_stats['events'])
        outf.write("Min: ", flow_stats['min'])
        outf.write("Max: ", flow_stats['max'])
        outf.write("Columns: ", flow_stats['columns'])
        outf.write("Markers: ", flow_stats['markers'])
        outf.write("Population: ", flow_stats['population'])
        outf.write("Population Freq: ", flow_stats['population_freq'])
        outf.write("Population Sample: ", flow_stats['population_sample'])
        outf.write("Population Per: ", flow_stats['population_per'])
        outf.write("Sample Data contains ", len(flow_stats['sample']), " events")
        outf.write("MIF_POP ", flow_stats['mfi_pop'])
