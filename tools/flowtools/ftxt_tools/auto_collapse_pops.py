#!/usr/bin/env python

######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################
from __future__ import print_function
import sys
import pandas as pd
from argparse import ArgumentParser


def auto_collapse(input_file, profile_file, output, report):
    profile_pop_list = {}
    pop_to_collapse = []
    markers = []
    with open(profile_file, "r") as pf:
        pffl = pf.readline()
        markers = pffl.strip().split("\t")
        for pfline in pf:
            line = pfline.strip().split("\t")
            pop = line[0]
            profil = "\t".join(line[1:-2])
            if profil in profile_pop_list:
                profile_pop_list[profil].append(pop)
            else:
                profile_pop_list[profil] = [pop]
    i = 1
    with open(report, "w") as rt:
        rt.write("New_Population\tFormer_Populations\t")
        rt.write("\t".join(markers[1:-2]) + "\n")
        for profs in profile_pop_list:
            pop_to_collapse.append(profile_pop_list[profs])
            pop_ls = ", ".join(profile_pop_list[profs])
            rt.write("\t".join([str(i), pop_ls, profs]) + "\n")
            i += 1
    df = pd.read_table(input_file, dtype={'Population': object})
    df['new_population'] = df.Population
    for i, sets_pop in enumerate(pop_to_collapse):
        df.loc[df['Population'].isin(sets_pop), ['new_population']] = i + 1

    df.Population = df.new_population
    df.drop(['new_population'], inplace=True, axis=1)
    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    parser = ArgumentParser(
             prog="auto_pop_collapse_from_profile",
             description="collapses FLOCK populations based on profile.")

    parser.add_argument(
            '-i',
            dest="input_file",
            required=True,
            help="FLOCK output file")

    parser.add_argument(
            '-o',
            dest="output",
            required=True,
            help="Name of the output file.")

    parser.add_argument(
            '-r',
            dest="report",
            required=True,
            help="Name of the report file.")

    parser.add_argument(
            '-p',
            dest="profile",
            required=True,
            help="File location for the profile.txt from FLOCK.")

    args = parser.parse_args()

    auto_collapse(args.input_file, args.profile, args.output, args.report)
