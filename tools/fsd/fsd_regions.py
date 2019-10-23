#!/usr/bin/env python

# Family size distribution of tags which were aligned to the reference genome
#
# Author: Monika Heinzl & Gundula Povysil, Johannes-Kepler University Linz (Austria)
# Contact: monika.heinzl@edumail.at
#
# Takes at least one TABULAR file with tags before the alignment to the SSCS,
# a BAM file with tags of reads that overlap the regions of the reference genome and
# an optional BED file with chromosome, start and stop position of the regions as input.
# The program produces a plot which shows the distribution of family sizes of the tags from the input files and
# a tabular file with the data of the plot.

# USAGE: python FSD_regions.py --inputFile filenameSSCS --inputName1 filenameSSCS
# --bamFile DCSbamFile --rangesFile BEDfile --output_tabular outptufile_name_tabular
# --output_pdf outputfile_name_pdf

import argparse
import collections
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pysam
from matplotlib.backends.backend_pdf import PdfPages

plt.switch_backend('agg')


def readFileReferenceFree(file, delim):
    with open(file, 'r') as dest_f:
        data_array = np.genfromtxt(dest_f, skip_header=0, delimiter=delim, comments='#', dtype='string')
        return data_array


def make_argparser():
    parser = argparse.ArgumentParser(description='Family Size Distribution of tags which were aligned to regions of the reference genome')
    parser.add_argument('--inputFile', help='Tabular File with three columns: ab or ba, tag and family size.')
    parser.add_argument('--inputName1')
    parser.add_argument('--bamFile', help='BAM file with aligned reads.')
    parser.add_argument('--rangesFile', default=None, help='BED file with chromosome, start and stop positions.')
    parser.add_argument('--output_pdf', default="data.pdf", type=str, help='Name of the pdf and tabular file.')
    parser.add_argument('--output_tabular', default="data.tabular", type=str, help='Name of the pdf and tabular file.')
    return parser


def compare_read_families_refGenome(argv):
    parser = make_argparser()
    args = parser.parse_args(argv[1:])

    firstFile = args.inputFile
    name1 = args.inputName1
    name1 = name1.split(".tabular")[0]
    bamFile = args.bamFile

    rangesFile = args.rangesFile
    title_file = args.output_pdf
    title_file2 = args.output_tabular
    sep = "\t"

    with open(title_file2, "w") as output_file, PdfPages(title_file) as pdf:
        data_array = readFileReferenceFree(firstFile, "\t")
        pysam.index(bamFile)

        bam = pysam.AlignmentFile(bamFile, "rb")
        qname_dict = collections.OrderedDict()

        if rangesFile != str(None):
            with open(rangesFile, 'r') as regs:
                range_array = np.genfromtxt(regs, skip_header=0, delimiter='\t', comments='#', dtype='string')

            if range_array.ndim == 0:
                print("Error: file has 0 lines")
                exit(2)

            if range_array.ndim == 1:
                chrList = range_array[0]
                start_posList = range_array[1].astype(int)
                stop_posList = range_array[2].astype(int)
                chrList = [chrList.tolist()]
                start_posList = [start_posList.tolist()]
                stop_posList = [stop_posList.tolist()]
            else:
                chrList = range_array[:, 0]
                start_posList = range_array[:, 1].astype(int)
                stop_posList = range_array[:, 2].astype(int)

            if len(start_posList) != len(stop_posList):
                print("start_positions and end_positions do not have the same length")
                exit(3)

            chrList = np.array(chrList)
            start_posList = np.array(start_posList).astype(int)
            stop_posList = np.array(stop_posList).astype(int)
            for chr, start_pos, stop_pos in zip(chrList, start_posList, stop_posList):
                chr_start_stop = "{}_{}_{}".format(chr, start_pos, stop_pos)
                qname_dict[chr_start_stop] = []
                for read in bam.fetch(chr.tobytes(), start_pos, stop_pos):
                    if not read.is_unmapped:
                        if re.search('_', read.query_name):
                            tags = re.split('_', read.query_name)[0]
                        else:
                            tags = read.query_name
                        qname_dict[chr_start_stop].append(tags)

        else:
            for read in bam.fetch():
                if not read.is_unmapped:
                    if re.search(r'_', read.query_name):
                        tags = re.split('_', read.query_name)[0]
                    else:
                        tags = read.query_name

                    if read.reference_name not in qname_dict.keys():
                        qname_dict[read.reference_name] = [tags]
                    else:
                        qname_dict[read.reference_name].append(tags)

        seq = np.array(data_array[:, 1])
        tags = np.array(data_array[:, 2])
        quant = np.array(data_array[:, 0]).astype(int)
        group = np.array(qname_dict.keys())

        all_ab = seq[np.where(tags == "ab")[0]]
        all_ba = seq[np.where(tags == "ba")[0]]
        quant_ab = quant[np.where(tags == "ab")[0]]
        quant_ba = quant[np.where(tags == "ba")[0]]

        seqDic_ab = dict(zip(all_ab, quant_ab))
        seqDic_ba = dict(zip(all_ba, quant_ba))

        lst_ab = []
        lst_ba = []
        quantAfterRegion = []
        length_regions = 0
        for i in group:
            lst_ab_r = []
            lst_ba_r = []
            seq_mut = qname_dict[i]
            if rangesFile == str(None):
                seq_mut, seqMut_index = np.unique(np.array(seq_mut), return_index=True)
            length_regions = length_regions + len(seq_mut) * 2
            for r in seq_mut:
                count_ab = seqDic_ab.get(r)
                count_ba = seqDic_ba.get(r)
                lst_ab_r.append(count_ab)
                lst_ab.append(count_ab)
                lst_ba_r.append(count_ba)
                lst_ba.append(count_ba)

            dataAB = np.array(lst_ab_r)
            dataBA = np.array(lst_ba_r)
            bigFamilies = np.where(dataAB > 20)[0]
            dataAB[bigFamilies] = 22
            bigFamilies = np.where(dataBA > 20)[0]
            dataBA[bigFamilies] = 22

            quantAll = np.concatenate((dataAB, dataBA))
            quantAfterRegion.append(quantAll)

        quant_ab = np.array(lst_ab)
        quant_ba = np.array(lst_ba)

        maximumX = np.amax(np.concatenate(quantAfterRegion))
        minimumX = np.amin(np.concatenate(quantAfterRegion))

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
        for i in range(0, len(group)):
            col.append(colors[i])

        counts = plt.hist(quantAfterRegion, bins=range(minimumX, maximumX + 1), stacked=False, label=group,
                          align="left", alpha=1, color=col, edgecolor="black", linewidth=1)
        ticks = np.arange(minimumX - 1, maximumX, 1)

        ticks1 = map(str, ticks)
        ticks1[len(ticks1) - 1] = ">20"
        plt.xticks(np.array(ticks), ticks1)
        count = np.bincount(map(int, quant_ab))  # original counts

        legend = "max. family size:\nabsolute frequency:\nrelative frequency:\n\ntotal nr. of reads:\n(before SSCS building)"
        plt.text(0.15, 0.085, legend, size=11, transform=plt.gcf().transFigure)

        legend = "AB\n{}\n{}\n{:.5f}\n\n{:,}".format(max(map(int, quant_ab)), count[len(count) - 1], float(count[len(count) - 1]) / sum(count), sum(np.array(data_array[:, 0]).astype(int)))
        plt.text(0.35, 0.105, legend, size=11, transform=plt.gcf().transFigure)

        count2 = np.bincount(map(int, quant_ba))  # original counts

        legend = "BA\n{}\n{}\n{:.5f}" \
            .format(max(map(int, quant_ba)), count2[len(count2) - 1], float(count2[len(count2) - 1]) / sum(count2))
        plt.text(0.45, 0.1475, legend, size=11, transform=plt.gcf().transFigure)

        plt.text(0.55, 0.2125, "total nr. of tags:", size=11, transform=plt.gcf().transFigure)
        plt.text(0.8, 0.2125, "{:,} ({:,})".format(length_regions, length_regions / 2), size=11,
                 transform=plt.gcf().transFigure)

        legend4 = "* In the plot, both family sizes of the ab and ba strands were used.\nWhereas the total numbers indicate only the single count of the tags per region.\n"
        plt.text(0.1, 0.01, legend4, size=11, transform=plt.gcf().transFigure)

        space = 0
        for i, count in zip(group, quantAfterRegion):
            plt.text(0.55, 0.15 - space, "{}:\n".format(i), size=11, transform=plt.gcf().transFigure)
            plt.text(0.8, 0.15 - space, "{:,}\n".format(len(count) / 2), size=11, transform=plt.gcf().transFigure)
            space = space + 0.02

        plt.legend(loc='upper right', fontsize=14, bbox_to_anchor=(0.9, 1), frameon=True)
        plt.xlabel("Family size", fontsize=14)
        plt.ylabel("Absolute Frequency", fontsize=14)
        plt.grid(b=True, which="major", color="#424242", linestyle=":")
        plt.margins(0.01, None)

        pdf.savefig(fig, bbox_inch="tight")
        plt.close()

        output_file.write("Dataset:{}{}\n".format(sep, name1))
        output_file.write("{}AB{}BA\n".format(sep, sep))
        output_file.write("max. family size:{}{}{}{}\n".format(sep, max(map(int, quant_ab)), sep, max(map(int, quant_ba))))
        output_file.write("absolute frequency:{}{}{}{}\n".format(sep, count[len(count) - 1], sep, count2[len(count2) - 1]))
        output_file.write("relative frequency:{}{:.3f}{}{:.3f}\n\n".format(sep, float(count[len(count) - 1]) / sum(count), sep, float(count2[len(count2) - 1]) / sum(count2)))
        output_file.write("total nr. of reads{}{}\n".format(sep, sum(np.array(data_array[:, 0]).astype(int))))
        output_file.write("total nr. of tags{}{} ({})\n".format(sep, length_regions, length_regions / 2))
        output_file.write("\n\nValues from family size distribution\n")
        output_file.write("{}".format(sep))
        for i in group:
            output_file.write("{}{}".format(i, sep))
        output_file.write("\n")

        j = 0
        for fs in counts[1][0:len(counts[1]) - 1]:
            if fs == 21:
                fs = ">20"
            else:
                fs = "={}".format(fs)
            output_file.write("FS{}{}".format(fs, sep))
            if len(group) == 1:
                output_file.write("{}{}".format(int(counts[0][j]), sep))
            else:
                for n in range(len(group)):
                    output_file.write("{}{}".format(int(counts[0][n][j]), sep))
            output_file.write("\n")
            j += 1
        output_file.write("sum{}".format(sep))

        if len(group) == 1:
            output_file.write("{}{}".format(int(sum(counts[0])), sep))
        else:
            for i in counts[0]:
                output_file.write("{}{}".format(int(sum(i)), sep))
        output_file.write("\n")
        output_file.write("\n\nIn the plot, both family sizes of the ab and ba strands were used.\nWhereas the total numbers indicate only the single count of the tags per region.\n")
        output_file.write("Region{}total nr. of tags per region\n".format(sep))
        for i, count in zip(group, quantAfterRegion):
            output_file.write("{}{}{}\n".format(i, sep, len(count) / 2))

    print("Files successfully created!")


if __name__ == '__main__':
    sys.exit(compare_read_families_refGenome(sys.argv))
