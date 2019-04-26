#!/usr/bin/env python

######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################
from __future__ import print_function
import sys
import os
from argparse import ArgumentParser
from jinja2 import Environment, FileSystemLoader

profile_key = {
    "1": "-",
    "2": "lo",
    "3": "+",
    "4": "hi"
}


def run_flowCL(phenotype, output_txt, output_pdf, tool):
    run_command = " ". join(["Rscript --slave --vanilla", tool, output_txt, phenotype])
    os.system(run_command)
    get_graph = " ".join(["mv flowCL_results/*.pdf", output_pdf])
    os.system(get_graph)
    return


def generate_flowCL_query(list_markers, list_types):
    if (len(list_markers) != len(list_types)):
        return("pb with headers")
    query = []
    # go through both lists, remove fsc/ssc
    for i in range(1, len(list_markers)):
        if not list_markers[i].startswith("FSC") and not list_markers[i].startswith("SSC"):
            query.append(list_markers[i].upper())
            query.append(profile_key[list_types[i]])
    # return concatenated string
    return("".join(query))


def translate_profiles(input_file, tool_dir, output, html_dir):
    os.mkdir(html_dir)

    tool = "/".join([tool_dir, "getOntology.R"])
    html_table = "".join([html_dir, "/CLprofiles.txt"])
    score_table = "".join(["cp ", input_file, " ", html_dir, "/scores.txt"])
    os.system(score_table)

    # read profile
    with open(input_file, "r") as flock_profiles, open(html_table, "w") as out:
        headers = flock_profiles.readline()
        headers = headers.strip()
        # get all headers except for last 2 (count + percentage)
        markers = headers.split("\t")[:-2]
        counter = 0

        out.write("Population\tFlowCL Query\tNb Results\tLink to PDF\t")
        out.write("Top Result Label\tTop Result Score\tTop Result CL\n")
        queries = {}
        # create marker query for each population
        for lines in flock_profiles:
            lines = lines.strip("\n")
            pop_profile = lines.split("\t")[:-2]
            flowcl_query = generate_flowCL_query(markers, pop_profile)
            counter += 1
            nb_results = "0"
            top_label = "no_match"
            top_score = "NA"
            top_CL = "NA"
            pdf_link = "NA"
            # check if query was run before
            if flowcl_query not in queries:
                # create filenames for results & graphs
                txt = "".join(["flowcl_pop", str(counter).zfill(2), ".txt"])
                text_result = "/".join([html_dir, txt])
                graph = "".join(["flowcl_pop", str(counter).zfill(2), ".pdf"])
                graph_output = "/".join([html_dir, graph])
                # run flowCL for each marker profile
                run_flowCL(flowcl_query, text_result, graph_output, tool)

                # test that text file exists if not results are all NAs:
                if os.path.isfile(text_result):
                    with open(text_result, "r") as res:
                        for line in res:
                            if line.startswith("Score"):
                                data = line.split(") ")
                                top_score = data[2][:-2]
                                tot_results = len(data) - 2
                                nb_results = str(tot_results)
                                if tot_results == 5:
                                    if len(data[6].split("+")) > 1:
                                        nb_results = "5+"
                            elif line.startswith("Cell ID"):
                                prep_link = line.split(") ")[1][:-2]
                                cl = prep_link.replace("_", ":")
                                link = "".join(['<a href="http://www.immport-labs.org/immport-ontology/public/home/home/', cl, '" target="_blank">'])
                                top_CL = "".join([link, prep_link, "</a>"])
                            elif line.startswith("Cell Label"):
                                top_label = line.split(") ")[1][:-2]
                                pdf_link = "".join(['<a href="', graph, '" target="_blank">PDF</a>'])
                                tmpflowcl_query = "".join(['<a href="', txt, '" target="_blank">', flowcl_query, '</a>'])

                    queries[flowcl_query] = {
                        "query": tmpflowcl_query,
                        "results": nb_results,
                        "pdf": pdf_link,
                        "label": top_label,
                        "score": top_score,
                        "CL": top_CL
                    }
            # write query results to CLprofiles.txt
            out.write("\t".join([pop_profile[0],
                                 queries[flowcl_query]["query"],
                                 queries[flowcl_query]["results"],
                                 queries[flowcl_query]["pdf"],
                                 queries[flowcl_query]["label"],
                                 queries[flowcl_query]["score"],
                                 queries[flowcl_query]["CL"]]) + "\n")

    env = Environment(loader=FileSystemLoader(tool_dir + "/templates"))
    template = env.get_template("profileCLs.template")

    real_directory = html_dir.replace("/job_working_directory", "")
    context = {'outputDirectory': real_directory}
    overview = template.render(**context)
    with open(output, "w") as outf:
        outf.write(overview)


if __name__ == "__main__":
    parser = ArgumentParser(
             prog="getCLs_from_profile",
             description="runs flowCL on a each population defined by FLOCK.")

    parser.add_argument(
            '-i',
            dest="input_file",
            required=True,
            help="File location for the profile.txt from FLOCK.")

    parser.add_argument(
            '-o',
            dest="output",
            required=True,
            help="Name of the output html file.")

    parser.add_argument(
            '-d',
            dest="html_dir",
            required=True,
            help="Path to html supporting directory.")

    parser.add_argument(
            '-t',
            dest="tool_dir",
            required=True,
            help="Path to the tool directory")

    args = parser.parse_args()

    translate_profiles(args.input_file, args.tool_dir, args.output, args.html_dir)
