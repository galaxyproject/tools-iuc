#!/usr/bin/env python
######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################
from __future__ import print_function
import sys
import os

from collections import defaultdict
from argparse import ArgumentParser
from jinja2 import Environment, FileSystemLoader


def generate_flowCL_query(list_markers, list_types):
    if (len(list_markers) != len(list_types)):
        return("pb with headers")
    query = []
    # go through both lists, remove fsc/ssc
    for i in range(0, len(list_markers)):
        if not list_markers[i].startswith("FSC") and not list_markers[i].startswith("SSC"):
            query.append(list_markers[i].upper())
            query.append(list_types[i])
    # return concatenated string
    return("".join(query))


def run_flowCL(phenotype, output_file, output_dir, tool_dir):
    os.mkdir(output_dir)
    tool = "/".join([tool_dir, "getOntology.R"])
    output_txt = "".join([output_dir, "/flowCL_run_summary.txt"])
    output_table = "".join([output_dir, "/flowCL_table.txt"])
    output_pdf = "".join([output_dir, "/flowCL_res.pdf"])
    run_command = " ". join(["Rscript --slave --vanilla", tool, output_txt, phenotype])
    os.system(run_command)

    table = defaultdict(list)
    labels = []
    nb_match = 0
    if os.path.isfile(output_txt):
        with open(output_txt, "r") as txt:
            check = txt.readline().strip()
            if (not check):
                sys.exit(2)
            else:
                i = -1
                for lines in txt:
                    data = lines.strip("\n").split("\"")
                    if data[0].strip():
                        labels.append(data[0].strip())
                        i += 1
                        if data[0].startswith("Score"):
                            count_matches = data[1].split(") ")
                            nb_match = len(count_matches) - 1
                    table[i].append(data[1])
    else:
        sys.stderr.write("There are no results with this query. Please check your markers if you believe there should be.")
        sys.exit(2)

    with open(output_table, "w") as tbl:
        tbl.write("1\t2\nQuery\t" + phenotype + "\n")
        for j in table:
            newline = " ".join(table[j])
            for k in range(1, nb_match + 1):
                cur_stg = "".join([str(k+1), ")"])
                new_stg = "".join(["<br>", cur_stg])
                newline = newline.replace(cur_stg, new_stg)

            if labels[j] == "Cell ID":
                cls = newline.split(" ")
                for m in range(0, len(cls)):
                    if cls[m].startswith("CL"):
                        cl_id = cls[m].replace("_", ":")
                        link = "".join(['<a href="http://www.immport-labs.org/immport-ontology/public/home/home/', cl_id, '" target="_blank">'])
                        cls[m] = "".join([link, cls[m], "</a>"])
                newline = " ".join(cls)
            tbl.write("\t".join([labels[j], newline]) + "\n")

    get_graph = " ".join(["mv flowCL_results/*.pdf", output_pdf])
    os.system(get_graph)

    env = Environment(loader=FileSystemLoader(tool_dir + "/templates"))
    template = env.get_template("flowCL.template")

    real_directory = output_dir.replace("/job_working_directory", "")
    context = {'outputDirectory': real_directory}
    overview = template.render(**context)
    with open(output_file, "w") as outf:
        outf.write(overview)
    return


if __name__ == "__main__":
    parser = ArgumentParser(
             prog="getOntology",
             description="runs flowCL on a set of markers.")

    parser.add_argument(
            '-m',
            dest="markers",
            required=True,
            action='append',
            help="marker queries.")

    parser.add_argument(
            '-y',
            dest="marker_types",
            required=True,
            action='append',
            help="marker queries.")

    parser.add_argument(
            '-o',
            dest="output_file",
            required=True,
            help="Name of the output html file.")

    parser.add_argument(
            '-d',
            dest="output_dir",
            required=True,
            help="Path to the html supporting directory")

    parser.add_argument(
            '-t',
            dest="tool_dir",
            required=True,
            help="Path to the tool directory")

    args = parser.parse_args()

    markers = [m.strip() for m in args.markers]
    query = generate_flowCL_query(markers, args.marker_types)
    run_flowCL(query, args.output_file, args.output_dir, args.tool_dir)
