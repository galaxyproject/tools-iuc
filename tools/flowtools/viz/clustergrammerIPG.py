#!/usr/bin/env python

######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################

# version 1.1 -- August 2017
# added checks for consistency between input files
# and upper limit on nb of cluster to look at
# version 1.2 -- June 2018
# added clustergrammer clustering parameters and normalization options


from __future__ import print_function
import sys
import os

from jinja2 import Environment, FileSystemLoader
import pandas as pd
from clustergrammer import Network


def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def get_indices(input_text):
    output_list = []
    default_value = ['i.e.:1,2,5', 'i.e.:3,6,7']
    if input_text in default_value:
        sys.exit(2)
    else:
        tmp_col = input_text.split(",")
        if len(tmp_col) == 1:
            if not tmp_col[0].strip():
                sys.exit(2)
            elif not is_integer(tmp_col[0].strip()):
                sys.exit(3)
            else:
                output_list.append(int(tmp_col[0].strip()) - 1)
        else:
            for c in range(0, len(tmp_col)):
                if not is_integer(tmp_col[c].strip()):
                    sys.exit(3)
                else:
                    output_list.append(int(tmp_col[c].strip()) - 1)
    return(output_list)


def prepare_heatmap(matrix_input, html_file, html_dir, tools_dir, categories, distance, linkage):
    # prepare directory and html
    os.mkdir(html_dir)

    env = Environment(loader=FileSystemLoader(tools_dir + "/templates"))
    template = env.get_template("clustergrammer.template")
    overview = template.render()
    with open(html_file, "w") as outf:
        outf.write(overview)

    json_output = html_dir + "/mult_view.json"

    net = Network()
    net.load_file(matrix_input)
    if (categories['row']):
        net.add_cats('row', categories['row'])
    if (categories['col']):
        net.add_cats('col', categories['col'])
    net.cluster(dist_type=distance, linkage_type=linkage)
    net.write_json_to_file('viz', json_output)


if __name__ == "__main__":

    args = sys.argv
    categories = {
        'row': [],
        'col': []
    }
    norm = {}

    if (len(args) > 7):
        df = pd.read_table(args[1])
        names = {
            'row': df.iloc[:, 0],
            'col': df.columns[1:]
        }

        tmp_string = "-=-".join(args[8:])
        print (tmp_string + "\n")
        # get categories
        cats = tmp_string.split("-=-new_cat-=-")
        for cat in cats:
            tmp_cat = cat.split("-=-")
            group = {
                "title": tmp_cat[1],
                "cats": {}
            }
            stg_groups = "--".join(tmp_cat[2:])
            groups = stg_groups.split("--new_label--")
            cat_indices = []
            for g in groups:
                print(g + "\n")
                elem = g.split("--")
                index_list = get_indices(elem[1])
                index_names = []
                for i in index_list:
                    if i in cat_indices:
                        sys.exit(4)
                    index_names.append(str(names[tmp_cat[0]][i]))
                cat_indices = cat_indices + index_list
                print(index_names, elem[0], sep="\t")
                group["cats"][elem[0]] = index_names
            categories[tmp_cat[0]].append(group)
            print(categories)


    prepare_heatmap(args[1], args[2], args[3], args[4], categories, args[5], args[6])
