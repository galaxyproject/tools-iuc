#!/usr/bin/env python
######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################

# version 1.1 - August 2017
# added upper limit to nb of clusters (40)
#
from __future__ import print_function
import sys
import os
import pandas as pd
import logging
import fileinput

from argparse import ArgumentParser
from jinja2 import Environment, FileSystemLoader
from collections import defaultdict

from color_palette import color_palette
from flowstatlib import gen_overview_stats
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


profile_key = {
    "1": "-",
    "2": "lo",
    "3": "+",
    "4": "hi"
}


# flow CL functions
def run_flowCL(phenotype, output_txt, output_pdf, tool):
    run_command = " ". join(["Rscript --slave --vanilla", tool, output_txt,
                             phenotype])
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


def translate_profiles(input_file, tool_dir, html_dir):
    tool = "/".join([tool_dir, "getOntology.R"])
    html_table = "".join([html_dir, "/CLprofiles.txt"])
    score_table = "".join(["cp ", input_file, " ", html_dir, "/scores.txt"])
    os.system(score_table)

    # read profile
    with open(input_file, "r") as flock_profiles, open(html_table, "w") as out:
        headers = flock_profiles.readline()
        headers = headers.strip()
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


# boxplots data massaging
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
                if (len(list(out[population][marker]))> 100):
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

# html generation
def gen_flow_overview(args):
    flow_stats = gen_overview_stats(args.input_file)
    if len(set(flow_stats['population'])) > 40:
        nbpop = str(len(set(flow_stats['population'])))
        sys.stderr.write("There are " + nbpop + " in the input file.")
        sys.exit(3)

    os.mkdir(args.output_directory)
    html_template = "genOverview.template"

    if args.scores:
        translate_profiles(args.scores, args.tool_directory, args.output_directory)
        html_template = "genOverviewCL.template"

    env = Environment(loader=FileSystemLoader(args.tool_directory + "/templates"))
    template = env.get_template(html_template)

    real_directory = args.output_directory.replace("/job_working_directory", "")
    context = {'outputDirectory': real_directory}
    overview = template.render(**context)
    with open(args.output_file, "w") as f:
        f.write(overview)

    flow_sample_file_name = args.output_directory + "/flow.sample"
    with open(flow_sample_file_name, "w") as flow_sample_file:
        flow_stats['sample'].to_csv(flow_sample_file, sep="\t", index=False, float_format='%.0f')

    flow_mfi_file_name = args.output_directory + "/flow.mfi"
    with open(flow_mfi_file_name, "w") as flow_mfi_file:
        flow_stats[args.mfi_calc].to_csv(flow_mfi_file, sep="\t", float_format='%.0f')

    mpop = "_".join([args.mfi_calc, "pop"])
    flow_mfi_pop_file_name = args.output_directory + "/flow.mfi_pop"
    with open(flow_mfi_pop_file_name, "w") as flow_mfi_pop_file:
        flow_stats[mpop].to_csv(flow_mfi_pop_file, sep="\t", index=False, float_format="%.2f")

    # box plot data
    boxplot_data = args.output_directory + "/boxplotData.json"
    resampled = get_boxplot_stats(args.input_file, flow_mfi_file_name, boxplot_data)

    # Generate the Images  -- eventually we should change that over to D3
    fcm = flow_stats['sample_data'].values
    colors = []
    for i, j in enumerate(flow_stats['sample_population']):
        colors.append(color_palette[j])

    for i in range(flow_stats['columns']):
        for j in range(flow_stats['columns']):
            file_name = "m" + str(i) + "_m" + str(j)
            ax = plt.subplot(1, 1, 1)
            plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.0, wspace=0.0, hspace=0.0)
            plt.scatter(fcm[:, i], fcm[:, j], s=1, c=colors, edgecolors='none')
            plt.axis([0, 1024, 0, 1024])
            plt.xticks([])
            plt.yticks([])
            F = plt.gcf()
            F.set_size_inches(1, 1)
            F.set_dpi(90)
            png_file = file_name + "_90X90.png"
            F.savefig(args.output_directory + "/" + png_file)
            plt.clf()

    flow_overview_file_name = args.output_directory + "/flow.overview"
    with open(flow_overview_file_name, "w") as flow_overview_file:
        flow_overview_file.write("<table>\n")
        flow_overview_file.write("<tr><td>&nbsp;</td>\n")
        for i in range(flow_stats['columns']):
            flow_overview_file.write("<td>" + flow_stats['markers'][i] + "</td>\n")

        for i in range(flow_stats['columns']):
            flow_overview_file.write("<tr>\n")
            flow_overview_file.write("<td>" + flow_stats['markers'][i] + "</td>\n")
            for j in range(flow_stats['columns']):
                file_name = "m" + str(j) + "_m" + str(i)
                image_file = file_name + "_90X90.png"
                flow_overview_file.write('<td><img src="' + image_file + '"/></td>')
            flow_overview_file.write("</tr>\n")

        flow_overview_file.write("</table>\n</body>\n<html>\n")

    if resampled:
        to_find = '<div id="outlierWarning" style="display:none;">'
        to_replace = '<div id="outlierWarning">'
        ## yay python 2.7
        ro = fileinput.input(args.output_file, inplace=True, backup=".bak")
        for roline in ro:
            print(roline.replace(to_find, to_replace), end='')
        ro.close()


if __name__ == "__main__":
    parser = ArgumentParser(
             prog="genOverview",
             description="Generate an overview plot of Flow results.")

    parser.add_argument(
            '-i',
            dest="input_file",
            required=True,
            help="File location for the Flow Text file.")

    parser.add_argument(
            '-o',
            dest="output_file",
            required=True,
            help="File location for the HTML output file.")

    parser.add_argument(
            '-d',
            dest="output_directory",
            required=True,
            help="Directory location for the Flow Plot.")

    parser.add_argument(
            '-M',
            dest="mfi_calc",
            required=True,
            help="what to calculate for centroids.")

    parser.add_argument(
            '-p',
            dest="scores",
            help="File location for FLOCK population scores.")

    parser.add_argument(
            '-t',
            dest="tool_directory",
            required=True,
            help="Location of the Tool Directory.")

    args = parser.parse_args()

    gen_flow_overview(args)
