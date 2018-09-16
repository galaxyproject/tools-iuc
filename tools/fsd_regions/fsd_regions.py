#!/usr/bin/env python

# Family size distribution of tags which were aligned to the reference genome
#
# Author: Monika Heinzl, Johannes-Kepler University Linz (Austria)
# Contact: monika.heinzl@edumail.at
#
# Takes at least one TABULAR file with tags before the alignment to the SSCS
# and a TXT with tags of reads that overlap the regions of the reference genome as input.
# The program produces a plot which shows the distribution of family sizes of the tags from the input files and
# a CSV file with the data of the plot.

# USAGE: python FSD_regions_1.6_FINAL.py --inputFile filenameSSCS --inputName1 filenameSSCS --ref_genome  filenameRefGenome --sep "characterWhichSeparatesCSVFile" --output_csv outptufile_name_csv --output_pdf outptufile_name_pdf

import argparse
import matplotlib.pyplot as plt
import numpy
import sys
from . import matplotlib
from .backends.backend_pdf import PdfPages

def readFileReferenceFree(file, delim):
    with open(file, 'r') as dest_f:
        data_array = numpy.genfromtxt(dest_f, skip_header=0, delimiter=delim, comments='#', dtype='string')
        return(data_array)

def make_argparser():
    parser = argparse.ArgumentParser(description='Family Size Distribution of tags which were aligned to regions of the reference genome')
    parser.add_argument('--inputFile',
                        help='Tabular File with three columns: ab or ba, tag and family size.')
    parser.add_argument('--inputName1')
    parser.add_argument('--ref_genome',
                        help='TXT File with tags of reads that overlap the region.')
    parser.add_argument('--output_pdf', default="data.pdf", type=str,
                       help='Name of the pdf and csv file.')
    parser.add_argument('--output_csv', default="data.csv", type=str,
                        help='Name of the pdf and csv file.')
    parser.add_argument('--sep', default=",",
                        help='Separator in the csv file.')
    return parser

def compare_read_families_refGenome(argv):
    parser = make_argparser()
    args = parser.parse_args(argv[1:])

    firstFile = args.inputFile
    name1 = args.inputName1
    name1 = name1.split(".tabular")[0]
    refGenome = args.ref_genome
    title_file = args.output_pdf
    title_file2 = args.output_csv
    sep = args.sep

    if type(sep) is not str or len(sep) > 1:
        print("Error: --sep must be a single character.")
        exit(3)

    with open(title_file2, "w") as output_file, PdfPages(title_file) as pdf:
        data_array = readFileReferenceFree(firstFile, "\t")

        mut_array = readFileReferenceFree(refGenome, " ")
        length_regions = len(mut_array)

        seq = numpy.array(data_array[:, 1])
        tags = numpy.array(data_array[:, 2])
        quant = numpy.array(data_array[:, 0]).astype(int)
        group = numpy.array(mut_array[:, 0])

        all_ab = seq[numpy.where(tags == "ab")[0]]
        all_ba = seq[numpy.where(tags == "ba")[0]]
        quant_ab = quant[numpy.where(tags == "ab")[0]]
        quant_ba = quant[numpy.where(tags == "ba")[0]]

        seqDic_ab = dict(zip(all_ab, quant_ab))
        seqDic_ba = dict(zip(all_ba, quant_ba))

        seq_mut = numpy.array(mut_array[:, 1])

        groupUnique, group_index = numpy.unique(group, return_index=True)
        groupUnique = groupUnique[numpy.argsort(group_index)]

        lst_ab = []
        for i in seq_mut:
            lst_ab.append(seqDic_ab.get(i))

        lst_ba = []
        for i in seq_mut:
            lst_ba.append(seqDic_ba.get(i))

        quant_ab = numpy.array(lst_ab)
        quant_ba = numpy.array(lst_ba)

        quantAfterRegion = []
        for i in groupUnique:
            dataAB = quant_ab[numpy.where(group == i)[0]]
            dataBA = quant_ba[numpy.where(group == i)[0]]
            bigFamilies = numpy.where(dataAB > 20)[0]
            dataAB[bigFamilies] = 22
            bigFamilies = numpy.where(dataBA > 20)[0]
            dataBA[bigFamilies] = 22

            quantAll = numpy.concatenate((dataAB, dataBA))
            quantAfterRegion.append(quantAll)

        maximumX = numpy.amax(numpy.concatenate(quantAfterRegion))
        minimumX = numpy.amin(numpy.concatenate(quantAfterRegion))

        # PLOT 
        plt.rc('figure', figsize=(11.69, 8.27))  # A4 format
        plt.rcParams['axes.facecolor'] = "E0E0E0"  # grey background color
        plt.rcParams['xtick.labelsize'] = 14
        plt.rcParams['ytick.labelsize'] = 14
        plt.rcParams['patch.edgecolor'] = "black"
        fig = plt.figure()
        plt.subplots_adjust(bottom=0.3)

        colors = ["#6E6E6E", "#0431B4", "#5FB404", "#B40431", "#F4FA58", "#DF7401", "#81DAF5"]

        col = []
        for i in range(0, len(groupUnique)):
            col.append(colors[i])

        counts = plt.hist(quantAfterRegion, bins=range(minimumX, maximumX + 1), stacked=False, label=groupUnique,
                          align="left", alpha=1, color=col, edgecolor="black", linewidth=1)
        ticks = numpy.arange(minimumX - 1, maximumX, 1)

        ticks1 = map(str, ticks)
        ticks1[len(ticks1) - 1] = ">20"
        plt.xticks(numpy.array(ticks), ticks1)
        count = numpy.bincount(map(int, quant_ab))  # original counts

        legend = "max. family size =\nabsolute frequency=\nrelative frequency=\n\ntotal nr. of reads="
        plt.text(0.15, 0.105, legend, size=11, transform=plt.gcf().transFigure)

        legend = "AB\n{}\n{}\n{:.5f}\n\n{:,}" \
            .format(max(map(int, quant_ab)), count[len(count) - 1], float(count[len(count) - 1]) / sum(count),
                    sum(numpy.array(data_array[:, 0]).astype(int)))
        plt.text(0.35, 0.105, legend, size=11, transform=plt.gcf().transFigure)

        count2 = numpy.bincount(map(int, quant_ba))  # original counts

        legend = "BA\n{}\n{}\n{:.5f}" \
            .format(max(map(int, quant_ba)), count2[len(count2) - 1], float(count2[len(count2) - 1]) / sum(count2))
        plt.text(0.45, 0.15, legend, size=11, transform=plt.gcf().transFigure)

        legend1 = "total nr. of tags="
        legend2 = "total numbers * \n{:,}".format(length_regions)
        plt.text(0.6, 0.2, legend1, size=11, transform=plt.gcf().transFigure)
        plt.text(0.75, 0.2, legend2, size=11, transform=plt.gcf().transFigure)
        legend4 = "* In the plot, both family sizes of the ab and ba strands were used.\nWhereas the total numbers indicate only the single count of the tags per region.\n"
        plt.text(0.1, 0.02, legend4, size=11, transform=plt.gcf().transFigure)

        space = numpy.arange(0, len(groupUnique), 0.02)
        for i, s, count in zip(groupUnique, space, quantAfterRegion):
            plt.text(0.6, 0.05 + s, "{}=\n".format(i), size=11, transform=plt.gcf().transFigure)
            plt.text(0.75, 0.05 + s, "{:,}\n".format(len(count) / 2), size=11, transform=plt.gcf().transFigure)

        plt.legend(loc='upper right', fontsize=14, bbox_to_anchor=(0.9, 1), frameon=True)
        # plt.title(name1, fontsize=14)
        plt.xlabel("Family size", fontsize=14)
        plt.ylabel("Absolute Frequency", fontsize=14)
        plt.grid(b=True, which="major", color="#424242", linestyle=":")
        plt.margins(0.01, None)

        # plt.savefig("{}_regions.pdf".format(title_file), bbox_inch="tight")
        pdf.savefig(fig, bbox_inch="tight")
        plt.close()

        output_file.write("Dataset:{}{}\n".format(sep, name1))
        output_file.write("{}AB{}BA\n".format(sep, sep))
        output_file.write("max. family size:{}{}{}{}\n".format(sep, max(map(int, quant_ab)), sep, max(map(int, quant_ba))))
        output_file.write("absolute frequency:{}{}{}{}\n".format(sep, count[len(count) - 1], sep, count2[len(count2) - 1]))
        output_file.write("relative frequency:{}{:.3f}{}{:.3f}\n\n".format(sep, float(count[len(count) - 1]) / sum(count), sep, float(count2[len(count2) - 1]) / sum(count2)))
        output_file.write("total nr. of reads{}{}\n".format(sep, sum(numpy.array(data_array[:, 0]).astype(int))))
        output_file.write("\n\nValues from family size distribution\n")
        output_file.write("{}".format(sep))
        for i in groupUnique:
            output_file.write("{}{}".format(i, sep))
        output_file.write("\n")
        j = 0
        for fs in counts[1][0:len(counts[1]) - 1]:
            if fs == 21:
                fs = ">20"
            else:
                fs = "={}".format(fs)
            output_file.write("FS{}{}".format(fs, sep))
            for n in range(len(groupUnique)):
                output_file.write("{}{}".format(int(counts[0][n][j]), sep))
            output_file.write("\n")
            j += 1
        output_file.write("sum{}".format(sep))
        for i in counts[0]:
            output_file.write("{}{}".format(int(sum(i)), sep))
        output_file.write("\n")
        output_file.write("\n\nIn the plot, both family sizes of the ab and ba strands were used.\nWhereas the total numbers indicate only the single count of the tags per region.\n")
        output_file.write("Region{}total nr. of tags per region\n".format(sep))
        for i, count in zip(groupUnique, quantAfterRegion):
            output_file.write("{}{}{}\n".format(i, sep, len(count) / 2))
        output_file.write("sum of tags{}{}\n".format(sep, length_regions))

    print("Files successfully created!")
    # print("Files saved under {}.pdf and {}.csv in {}!".format(title_file, title_file, os.getcwd()))


if __name__ == '__main__':
    sys.exit(compare_read_families_refGenome(sys.argv))
