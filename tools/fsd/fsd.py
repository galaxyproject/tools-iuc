#!/usr/bin/env python

# Family size distribution of SSCSs
#
# Author: Monika Heinzl, Johannes-Kepler University Linz (Austria)
# Contact: monika.heinzl@edumail.at
#
# Takes at least one TABULAR file with tags before the alignment to the SSCS, but up to 4 files can be provided, as input.
# The program produces a plot which shows the distribution of family sizes of the all SSCSs from the input files and
# a tabular file with the data of the plot, as well as a TXT file with all tags of the DCS and their family sizes.
# If only one file is provided, then a family size distribution, which is separated after SSCSs without a partner and DCSs, is produced.
# Whereas a family size distribution with multiple data in one plot is produced, when more than one file (up to 4) is given.

# USAGE: python FSD_Galaxy_1.4_commandLine_FINAL.py --inputFile1 filename --inputName1 filename --inputFile2 filename2 --inputName2 filename2 --inputFile3 filename3 --inputName3 filename3 --inputFile4 filename4 --inputName4 filename4 --log_axis --output_tabular outptufile_name_tabular --output_pdf outptufile_name_pdf

import argparse
import sys

import matplotlib.pyplot as plt
import numpy
from matplotlib.backends.backend_pdf import PdfPages

plt.switch_backend('agg')


def readFileReferenceFree(file):
    with open(file, 'r') as dest_f:
        data_array = numpy.genfromtxt(dest_f, skip_header=0, delimiter='\t', comments='#', dtype='string')
        return(data_array)


def make_argparser():
    parser = argparse.ArgumentParser(description='Family Size Distribution of duplex sequencing data')
    parser.add_argument('--inputFile1', help='Tabular File with three columns: ab or ba, tag and family size.')
    parser.add_argument('--inputName1')
    parser.add_argument('--inputFile2', default=None, help='Tabular File with three columns: ab or ba, tag and family size.')
    parser.add_argument('--inputName2')
    parser.add_argument('--inputFile3', default=None, help='Tabular File with three columns: ab or ba, tag and family size.')
    parser.add_argument('--inputName3')
    parser.add_argument('--inputFile4', default=None, help='Tabular File with three columns: ab or ba, tag and family size.')
    parser.add_argument('--inputName4')
    parser.add_argument('--log_axis', action="store_false", help='Transform y axis in log scale.')
    parser.add_argument('--rel_freq', action="store_false", help='If False, the relative frequencies are displayed.')
    parser.add_argument('--output_pdf', default="data.pdf", type=str, help='Name of the pdf file.')
    parser.add_argument('--output_tabular', default="data.tabular", type=str, help='Name of the tabular file.')
    return parser


def compare_read_families(argv):

    parser = make_argparser()
    args = parser.parse_args(argv[1:])

    firstFile = args.inputFile1
    name1 = args.inputName1

    secondFile = args.inputFile2
    name2 = args.inputName2
    thirdFile = args.inputFile3
    name3 = args.inputName3
    fourthFile = args.inputFile4
    name4 = args.inputName4
    log_axis = args.log_axis
    rel_freq = args.rel_freq

    title_file = args.output_tabular
    title_file2 = args.output_pdf

    sep = "\t"

    plt.rc('figure', figsize=(11.69, 8.27))  # A4 format
    plt.rcParams['patch.edgecolor'] = "black"
    plt.rcParams['axes.facecolor'] = "E0E0E0"  # grey background color
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    list_to_plot = []
    label = []
    data_array_list = []
    list_to_plot_original = []
    colors = []
    bins = numpy.arange(1, 22)
    with open(title_file, "w") as output_file, PdfPages(title_file2) as pdf:
        fig = plt.figure()
        fig.subplots_adjust(left=0.12, right=0.97, bottom=0.23, top=0.94, hspace=0)
        fig2 = plt.figure()
        fig2.subplots_adjust(left=0.12, right=0.97, bottom=0.23, top=0.94, hspace=0)

        # plt.subplots_adjust(bottom=0.25)
        if firstFile is not None:
            file1 = readFileReferenceFree(firstFile)
            integers = numpy.array(file1[:, 0]).astype(int)  # keep original family sizes
            list_to_plot_original.append(integers)
            colors.append("#0000FF")

            # for plot: replace all big family sizes by 22
            # data1 = numpy.array(file1[:, 0]).astype(int)
            # bigFamilies = numpy.where(data1 > 20)[0]
            # data1[bigFamilies] = 22
            data1 = numpy.clip(integers, bins[0], bins[-1])
            name1 = name1.split(".tabular")[0]
            if len(name1) > 40:
                name1 = name1[:40]
            list_to_plot.append(data1)
            label.append(name1)
            data_array_list.append(file1)

            legend = "\n\n\n{}".format(name1)
            fig.text(0.05, 0.11, legend, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.05, 0.11, legend, size=10, transform=plt.gcf().transFigure)

            legend1 = "singletons:\nnr. of tags\n{:,} ({:.3f})".format(numpy.bincount(data1)[1],
                                                                       float(numpy.bincount(data1)[1]) / len(data1))
            fig.text(0.32, 0.11, legend1, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.32, 0.11, legend1, size=10, transform=plt.gcf().transFigure)

            legend3b = "PE reads\n{:,} ({:.3f})".format(numpy.bincount(data1)[1],
                                                        float(numpy.bincount(data1)[1]) / sum(integers))
            fig.text(0.45, 0.11, legend3b, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.45, 0.11, legend3b, size=10, transform=plt.gcf().transFigure)

            legend4 = "family size > 20:\nnr. of tags\n{:,} ({:.3f})".format(len(integers[integers > 20]),
                                                                             float(len(integers[integers > 20])) / len(integers))
            fig.text(0.58, 0.11, legend4, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.58, 0.11, legend4, size=10, transform=plt.gcf().transFigure)

            legend5 = "PE reads\n{:,} ({:.3f})".format(sum(integers[integers > 20]),
                                                       float(sum(integers[integers > 20])) / sum(integers))
            fig.text(0.70, 0.11, legend5, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.70, 0.11, legend5, size=10, transform=plt.gcf().transFigure)

            legend6 = "total nr. of\ntags\n{:,}".format(len(data1))
            fig.text(0.82, 0.11, legend6, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.82, 0.11, legend6, size=10, transform=plt.gcf().transFigure)

            legend6b = "PE reads\n{:,}".format(sum(integers))
            fig.text(0.89, 0.11, legend6b, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.89, 0.11, legend6b, size=10, transform=plt.gcf().transFigure)

        if secondFile != str(None):
            file2 = readFileReferenceFree(secondFile)
            integers2 = numpy.array(file2[:, 0]).astype(int)  # keep original family sizes
            list_to_plot_original.append(integers2)
            colors.append("#298A08")

            # data2 = numpy.asarray(file2[:, 0]).astype(int)
            # bigFamilies2 = numpy.where(data2 > 20)[0]
            # data2[bigFamilies2] = 22

            data2 = numpy.clip(integers2, bins[0], bins[-1])
            list_to_plot.append(data2)
            name2 = name2.split(".tabular")[0]
            if len(name2) > 40:
                name2 = name2[:40]
            label.append(name2)
            data_array_list.append(file2)

            fig.text(0.05, 0.09, name2, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.05, 0.09, name2, size=10, transform=plt.gcf().transFigure)

            legend1 = "{:,} ({:.3f})".format(numpy.bincount(data2)[1], float(numpy.bincount(data2)[1]) / len(data2))
            fig.text(0.32, 0.09, legend1, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.32, 0.09, legend1, size=10, transform=plt.gcf().transFigure)

            legend3 = "{:,} ({:.3f})".format(numpy.bincount(data2)[1], float(numpy.bincount(data2)[1]) / sum(integers2))
            fig.text(0.45, 0.09, legend3, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.45, 0.09, legend3, size=10, transform=plt.gcf().transFigure)

            legend4 = "{:,} ({:.3f})".format(len(integers2[integers2 > 20]),
                                             float(len(integers2[integers2 > 20])) / len(integers2))
            fig.text(0.58, 0.09, legend4, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.58, 0.09, legend4, size=10, transform=plt.gcf().transFigure)

            legend5 = "{:,} ({:.3f})".format(sum(integers2[integers2 > 20]),
                                             float(sum(integers2[integers2 > 20])) / sum(integers2))
            fig.text(0.70, 0.09, legend5, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.70, 0.09, legend5, size=10, transform=plt.gcf().transFigure)

            legend6 = "{:,}".format(len(data2))
            fig.text(0.82, 0.09, legend6, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.82, 0.09, legend6, size=10, transform=plt.gcf().transFigure)

            legend6b = "{:,}".format(sum(integers2))
            fig.text(0.89, 0.09, legend6b, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.89, 0.09, legend6b, size=10, transform=plt.gcf().transFigure)

        if thirdFile is not None:
            file3 = readFileReferenceFree(thirdFile)
            integers3 = numpy.array(file3[:, 0]).astype(int)  # keep original family sizes
            list_to_plot_original.append(integers3)
            colors.append("#DF0101")

            # data3 = numpy.asarray(file3[:, 0]).astype(int)
            # bigFamilies3 = numpy.where(data3 > 20)[0]
            # data3[bigFamilies3] = 22

            data3 = numpy.clip(integers3, bins[0], bins[-1])
            list_to_plot.append(data3)
            name3 = name3.split(".tabular")[0]
            if len(name3) > 40:
                name3 = name3[:40]
            label.append(name3)
            data_array_list.append(file3)

            fig.text(0.05, 0.07, name3, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.05, 0.07, name3, size=10, transform=plt.gcf().transFigure)

            legend1 = "{:,} ({:.3f})".format(numpy.bincount(data3)[1], float(numpy.bincount(data3)[1]) / len(data3))
            fig.text(0.32, 0.07, legend1, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.32, 0.07, legend1, size=10, transform=plt.gcf().transFigure)

            legend3b = "{:,} ({:.3f})".format(numpy.bincount(data3)[1],
                                              float(numpy.bincount(data3)[1]) / sum(integers3))
            fig.text(0.45, 0.07, legend3b, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.45, 0.07, legend3b, size=10, transform=plt.gcf().transFigure)

            legend4 = "{:,} ({:.3f})".format(len(integers3[integers3 > 20]),
                                             float(len(integers3[integers3 > 20])) / len(integers3))
            fig.text(0.58, 0.07, legend4, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.58, 0.07, legend4, size=10, transform=plt.gcf().transFigure)

            legend5 = "{:,} ({:.3f})".format(sum(integers3[integers3 > 20]),
                                             float(sum(integers3[integers3 > 20])) / sum(integers3))
            fig.text(0.70, 0.07, legend5, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.70, 0.07, legend5, size=10, transform=plt.gcf().transFigure)

            legend6 = "{:,}".format(len(data3))
            fig.text(0.82, 0.07, legend6, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.82, 0.07, legend6, size=10, transform=plt.gcf().transFigure)

            legend6b = "{:,}".format(sum(integers3))
            fig.text(0.89, 0.07, legend6b, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.89, 0.07, legend6b, size=10, transform=plt.gcf().transFigure)

        if fourthFile is not None:
            file4 = readFileReferenceFree(fourthFile)
            integers4 = numpy.array(file4[:, 0]).astype(int)  # keep original family sizes
            list_to_plot_original.append(integers4)
            colors.append("#04cec7")

            # data4 = numpy.asarray(file4[:, 0]).astype(int)
            # bigFamilies4 = numpy.where(data4 > 20)[0]
            # data4[bigFamilies4] = 22
            data4 = numpy.clip(integers4, bins[0], bins[-1])
            list_to_plot.append(data4)
            name4 = name4.split(".tabular")[0]
            if len(name4) > 40:
                name4 = name4[:40]
            label.append(name4)
            data_array_list.append(file4)

            fig.text(0.05, 0.05, name4, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.05, 0.05, name4, size=10, transform=plt.gcf().transFigure)

            legend1 = "{:,} ({:.3f})".format(numpy.bincount(data4)[1], float(numpy.bincount(data4)[1]) / len(data4))
            fig.text(0.32, 0.05, legend1, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.32, 0.05, legend1, size=10, transform=plt.gcf().transFigure)

            legend3b = "{:,} ({:.3f})".format(numpy.bincount(data4)[1],
                                              float(numpy.bincount(data4)[1]) / sum(integers4))
            fig.text(0.45, 0.05, legend3b, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.45, 0.05, legend3b, size=10, transform=plt.gcf().transFigure)

            legend4 = "{:,} ({:.3f})".format(len(integers4[integers4 > 20]),
                                             float(len(integers4[integers4 > 20])) / len(integers4))
            fig.text(0.58, 0.05, legend4, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.58, 0.05, legend4, size=10, transform=plt.gcf().transFigure)

            legend5 = "{:,} ({:.3f})".format(sum(integers4[integers4 > 20]),
                                             float(sum(integers4[integers4 > 20])) / sum(integers4))
            fig.text(0.70, 0.05, legend5, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.70, 0.05, legend5, size=10, transform=plt.gcf().transFigure)

            legend6 = "{:,}".format(len(data4))
            fig.text(0.82, 0.05, legend6, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.82, 0.05, legend6, size=10, transform=plt.gcf().transFigure)

            legend6b = "{:,}".format(sum(integers4))
            fig.text(0.89, 0.05, legend6b, size=10, transform=plt.gcf().transFigure)
            fig2.text(0.89, 0.05, legend6b, size=10, transform=plt.gcf().transFigure)

        # maximumX = numpy.amax(numpy.concatenate(list_to_plot))
        # minimumX = numpy.amin(numpy.concatenate(list_to_plot))
        list_to_plot2 = list_to_plot

        if rel_freq:
            ylab = "Relative Frequency"
        else:
            ylab = "Absolute Frequency"

        # PLOT FSD based on tags
        fig.suptitle('Family Size Distribution (FSD) based on families', fontsize=14)
        ax = fig.add_subplot(1, 1, 1)
        ticks = numpy.arange(1, 22, 1)
        ticks1 = map(str, ticks)
        ticks1[len(ticks1) - 1] = ">20"
        ax.set_xticks([], [])
        if rel_freq:
            w = [numpy.zeros_like(data) + 1. / len(data) for data in list_to_plot2]
            counts = ax.hist(list_to_plot2, weights=w, bins=numpy.arange(1, 23), stacked=False, edgecolor="black", color=colors, linewidth=1, label=label, align="left", alpha=0.7, rwidth=0.8)
            ax.set_ylim(0, 1.07)
        else:
            counts = ax.hist(list_to_plot2, bins=numpy.arange(1, 23), stacked=False, edgecolor="black", linewidth=1, label=label, align="left", alpha=0.7, rwidth=0.8, color=colors)
        ax.set_xticks(numpy.array(ticks))
        ax.set_xticklabels(ticks1)
        ax.legend(loc='upper right', fontsize=14, frameon=True, bbox_to_anchor=(0.9, 1))

        ax.set_ylabel(ylab, fontsize=14)
        ax.set_xlabel("Family size", fontsize=14)
        if log_axis:
            ax.set_yscale('log')
        ax.grid(b=True, which="major", color="#424242", linestyle=":")
        ax.margins(0.01, None)
        pdf.savefig(fig)

        # PLOT FSD based on PE reads
        fig2.suptitle('Family Size Distribution (FSD) based on PE reads', fontsize=14)
        ax2 = fig2.add_subplot(1, 1, 1)
        ticks = numpy.arange(1, 22)
        ticks1 = map(str, ticks)
        ticks1[len(ticks1) - 1] = ">20"
        reads = []
        reads_rel = []

        barWidth = 0 - (len(list_to_plot) + 1) / 2 * 1. / (len(list_to_plot) + 1)
        ax2.set_xticks([], [])

        for i in range(len(list_to_plot2)):
            x = list(numpy.arange(1, 22).astype(float))
            unique, c = numpy.unique(list_to_plot2[i], return_counts=True)
            y = unique * c
            if sum(list_to_plot_original[i] > 20) > 0:
                y[len(y) - 1] = sum(list_to_plot_original[i][list_to_plot_original[i] > 20])
            y = [y[x[idx] == unique][0] if x[idx] in unique else 0 for idx in range(len(x))]
            reads.append(y)
            reads_rel.append(list(numpy.float_(y)) / sum(y))

            if len(list_to_plot2) == 1:
                x = [xi * 0.5 for xi in x]
                w = 0.4
            else:
                x = [xi + barWidth for xi in x]
                w = 1. / (len(list_to_plot) + 1)
            if rel_freq:
                ax2.bar(x, list(numpy.float_(y)) / numpy.sum(y), align="edge", width=w, edgecolor="black", label=label[i], linewidth=1, alpha=0.7, color=colors[i])
                ax2.set_ylim(0, 1.07)
            else:
                ax2.bar(x, y, align="edge", width=w, edgecolor="black", label=label[i], linewidth=1, alpha=0.7, color=colors[i])
            if i == len(list_to_plot2) - 1:
                barWidth += 1. / (len(list_to_plot) + 1) + 1. / (len(list_to_plot) + 1)
            else:
                barWidth += 1. / (len(list_to_plot) + 1)

        ax2.legend(loc='upper right', fontsize=14, frameon=True, bbox_to_anchor=(0.9, 1))

        if len(list_to_plot2) == 1:
            ax2.set_xticks(numpy.array([xi + 0.2 for xi in x]))
        else:
            ax2.set_xticks(numpy.array(ticks))
        ax2.set_xticklabels(ticks1)
        ax2.set_xlabel("Family size", fontsize=14)
        ax2.set_ylabel(ylab, fontsize=14)
        if log_axis:
            ax2.set_yscale('log')
        ax2.grid(b=True, which="major", color="#424242", linestyle=":")
        ax2.margins(0.01, None)

        pdf.savefig(fig2)
        plt.close()

        # write data to CSV file tags
        counts = [numpy.bincount(di, minlength=22)[1:] for di in list_to_plot2]  # original counts of family sizes
        output_file.write("Values from family size distribution with all datasets based on families\n")
        output_file.write("\nFamily size")
        for i in label:
            output_file.write("{}{}".format(sep, i))
        output_file.write("\n")
        j = 0
        for fs in bins:
            if fs == 21:
                fs = ">20"
            else:
                fs = "={}".format(fs)
            output_file.write("FS{}{}".format(fs, sep))
            for n in range(len(label)):
                output_file.write("{}{}".format(int(counts[n][j]), sep))
            output_file.write("\n")
            j += 1
        output_file.write("sum{}".format(sep))
        for i in counts:
            output_file.write("{}{}".format(int(sum(i)), sep))

        # write data to CSV file PE reads
        output_file.write("\n\nValues from family size distribution with all datasets based on PE reads\n")
        output_file.write("\nFamily size")
        for i in label:
            output_file.write("{}{}".format(sep, i))
        output_file.write("\n")
        j = 0

        for fs in bins:
            if fs == 21:
                fs = ">20"
            else:
                fs = "={}".format(fs)
            output_file.write("FS{}{}".format(fs, sep))
            if len(label) == 1:
                output_file.write("{}{}".format(int(reads[0][j]), sep))
            else:
                for n in range(len(label)):
                    output_file.write("{}{}".format(int(reads[n][j]), sep))
            output_file.write("\n")
            j += 1
        output_file.write("sum{}".format(sep))
        if len(label) == 1:
            output_file.write("{}{}".format(int(sum(numpy.concatenate(reads))), sep))
        else:
            for i in reads:
                output_file.write("{}{}".format(int(sum(i)), sep))
        output_file.write("\n")

        # Family size distribution after DCS and SSCS
        for dataset, data_o, name_file in zip(list_to_plot, data_array_list, label):
            # maximumX = numpy.amax(dataset)
            # minimumX = numpy.amin(dataset)

            tags = numpy.array(data_o[:, 2])
            seq = numpy.array(data_o[:, 1])
            data = numpy.array(dataset)
            data_o = numpy.array(data_o[:, 0]).astype(int)
            # find all unique tags and get the indices for ALL tags, but only once
            u, index_unique, c = numpy.unique(numpy.array(seq), return_counts=True, return_index=True)
            d = u[c > 1]

            # get family sizes, tag for duplicates
            duplTags_double = data[numpy.in1d(seq, d)]
            duplTags_double_o = data_o[numpy.in1d(seq, d)]

            duplTags = duplTags_double[0::2]  # ab of DCS
            duplTags_o = duplTags_double_o[0::2]  # ab of DCS

            duplTagsBA = duplTags_double[1::2]  # ba of DCS
            duplTagsBA_o = duplTags_double_o[1::2]  # ba of DCS

            # duplTags_double_tag = tags[numpy.in1d(seq, d)]
            # duplTags_double_seq = seq[numpy.in1d(seq, d)]

            # get family sizes for SSCS with no partner
            ab = numpy.where(tags == "ab")[0]
            abSeq = seq[ab]
            ab_o = data_o[ab]
            ab = data[ab]

            ba = numpy.where(tags == "ba")[0]
            baSeq = seq[ba]
            ba_o = data_o[ba]
            ba = data[ba]

            dataAB = ab[numpy.in1d(abSeq, d, invert=True)]
            dataAB_o = ab_o[numpy.in1d(abSeq, d, invert=True)]

            dataBA = ba[numpy.in1d(baSeq, d, invert=True)]
            dataBA_o = ba_o[numpy.in1d(baSeq, d, invert=True)]

            list1 = [duplTags_double, dataAB, dataBA]  # list for plotting
            list1_o = [duplTags_double_o, dataAB_o, dataBA_o]  # list for plotting

            # information for family size >= 3
            dataAB_FS3 = dataAB[dataAB >= 3]
            dataAB_FS3_o = dataAB_o[dataAB_o >= 3]
            dataBA_FS3 = dataBA[dataBA >= 3]
            dataBA_FS3_o = dataBA_o[dataBA_o >= 3]
            # ab_FS3 = ab[ab >= 3]
            # ba_FS3 = ba[ba >= 3]
            # ab_FS3_o = ab_o[ab_o >= 3]
            # ba_FS3_o = ba_o[ba_o >= 3]

            duplTags_FS3 = duplTags[(duplTags >= 3) & (duplTagsBA >= 3)]  # ab+ba with FS>=3
            duplTags_FS3_BA = duplTagsBA[(duplTags >= 3) & (duplTagsBA >= 3)]  # ba+ab with FS>=3
            duplTags_double_FS3 = len(duplTags_FS3) + len(duplTags_FS3_BA)  # both ab and ba strands with FS>=3

            # original FS
            duplTags_FS3_o = duplTags_o[(duplTags_o >= 3) & (duplTagsBA_o >= 3)]  # ab+ba with FS>=3
            duplTags_FS3_BA_o = duplTagsBA_o[(duplTags_o >= 3) & (duplTagsBA_o >= 3)]  # ba+ab with FS>=3
            duplTags_double_FS3_o = sum(duplTags_FS3_o) + sum(duplTags_FS3_BA_o)  # both ab and ba strands with FS>=3

            fig = plt.figure()
            plt.subplots_adjust(left=0.12, right=0.97, bottom=0.3, top=0.94, hspace=0)

            if rel_freq:
                w = [numpy.zeros_like(dj) + 1. / len(numpy.concatenate(list1)) for dj in list1]
                plt.hist(list1, bins=numpy.arange(1, 23), stacked=True, label=["duplex", "ab", "ba"], weights=w, edgecolor="black", linewidth=1, align="left", color=["#FF0000", "#5FB404", "#FFBF00"], rwidth=0.8)
                plt.ylim(0, 1.07)
            else:
                plt.hist(list1, bins=numpy.arange(1, 23), stacked=True, label=["duplex", "ab", "ba"], edgecolor="black", linewidth=1, align="left", color=["#FF0000", "#5FB404", "#FFBF00"], rwidth=0.8)

            # tick labels of x axis
            ticks = numpy.arange(1, 22, 1)
            ticks1 = map(str, ticks)
            ticks1[len(ticks1) - 1] = ">20"
            plt.xticks(numpy.array(ticks), ticks1)
            # singl = counts[0][2][0]  # singletons
            singl = len(data_o[data_o == 1])
            last = len(data_o[data_o > 20])  # large families
            if log_axis:
                plt.yscale('log')
            plt.legend(loc='upper right', fontsize=14, bbox_to_anchor=(0.9, 1), frameon=True)
            plt.title("{}: FSD based on families".format(name_file), fontsize=14)
            plt.xlabel("Family size", fontsize=14)
            plt.ylabel(ylab, fontsize=14)
            plt.margins(0.01, None)
            plt.grid(b=True, which="major", color="#424242", linestyle=":")

            # extra information beneath the plot
            legend = "SSCS ab= \nSSCS ba= \nDCS (total)= \ntotal nr. of tags="
            plt.text(0.1, 0.09, legend, size=10, transform=plt.gcf().transFigure)

            legend = "nr. of tags\n\n{:,}\n{:,}\n{:,} ({:,})\n{:,} ({:,})".format(len(dataAB), len(dataBA),
                                                                                  len(duplTags), len(duplTags_double), (len(dataAB) + len(dataBA) + len(duplTags)),
                                                                                  (len(ab) + len(ba)))
            plt.text(0.23, 0.09, legend, size=10, transform=plt.gcf().transFigure)

            legend5 = "PE reads\n\n{:,}\n{:,}\n{:,} ({:,})\n{:,} ({:,})".format(sum(dataAB_o), sum(dataBA_o),
                                                                                sum(duplTags_o), sum(duplTags_double_o),
                                                                                (sum(dataAB_o) + sum(dataBA_o) + sum(duplTags_o)),
                                                                                (sum(ab_o) + sum(ba_o)))
            plt.text(0.38, 0.09, legend5, size=10, transform=plt.gcf().transFigure)

            legend = "rel. freq. of tags\nunique\n{:.3f}\n{:.3f}\n{:.3f}\n{:,}".format(
                float(len(dataAB)) / (len(dataAB) + len(dataBA) + len(duplTags)),
                float(len(dataBA)) / (len(dataAB) + len(dataBA) + len(duplTags)),
                float(len(duplTags)) / (len(dataAB) + len(dataBA) + len(duplTags)),
                (len(dataAB) + len(dataBA) + len(duplTags)))
            plt.text(0.54, 0.09, legend, size=10, transform=plt.gcf().transFigure)

            legend = "total\n{:.3f}\n{:.3f}\n{:.3f} ({:.3f})\n{:,}".format(float(len(dataAB)) / (len(ab) + len(ba)),
                                                                           float(len(dataBA)) / (len(ab) + len(ba)),
                                                                           float(len(duplTags)) / (len(ab) + len(ba)),
                                                                           float(len(duplTags_double)) / (len(ab) + len(ba)),
                                                                           (len(ab) + len(ba)))
            plt.text(0.64, 0.09, legend, size=10, transform=plt.gcf().transFigure)

            legend1 = "\nsingletons:\nfamily size > 20:"
            plt.text(0.1, 0.03, legend1, size=10, transform=plt.gcf().transFigure)

            legend4 = "{:,}\n{:,}".format(singl, last)
            plt.text(0.23, 0.03, legend4, size=10, transform=plt.gcf().transFigure)
            legend3 = "{:.3f}\n{:.3f}".format(float(singl) / len(data), float(last) / len(data))
            plt.text(0.64, 0.03, legend3, size=10, transform=plt.gcf().transFigure)

            legend3 = "\n\n{:,}".format(sum(data_o[data_o > 20]))
            plt.text(0.38, 0.03, legend3, size=10, transform=plt.gcf().transFigure)

            legend3 = "{:.3f}\n{:.3f}".format(float(singl) / sum(data_o), float(sum(data_o[data_o > 20])) / sum(data_o))
            plt.text(0.84, 0.03, legend3, size=10, transform=plt.gcf().transFigure)

            legend = "PE reads\nunique\n{:.3f}\n{:.3f}\n{:.3f}\n{:,}".format(
                float(sum(dataAB_o)) / (sum(dataAB_o) + sum(dataBA_o) + sum(duplTags_o)),
                float(sum(dataBA_o)) / (sum(dataAB_o) + sum(dataBA_o) + sum(duplTags_o)),
                float(sum(duplTags_o)) / (sum(dataAB_o) + sum(dataBA_o) + sum(duplTags_o)),
                (sum(dataAB_o) + sum(dataBA_o) + sum(duplTags_o)))
            plt.text(0.74, 0.09, legend, size=10, transform=plt.gcf().transFigure)

            legend = "total\n{:.3f}\n{:.3f}\n{:.3f} ({:.3f})\n{:,}".format(
                float(sum(dataAB_o)) / (sum(ab_o) + sum(ba_o)),
                float(sum(dataBA_o)) / (sum(ab_o) + sum(ba_o)),
                float(sum(duplTags_o)) / (sum(ab_o) + sum(ba_o)),
                float(sum(duplTags_double_o)) / (sum(ab_o) + sum(ba_o)), (sum(ab_o) + sum(ba_o)))
            plt.text(0.84, 0.09, legend, size=10, transform=plt.gcf().transFigure)

            pdf.savefig(fig)
            plt.close()

            # PLOT FSD based on PE reads
            fig3 = plt.figure()
            plt.subplots_adjust(left=0.12, right=0.97, bottom=0.3, top=0.94, hspace=0)

            fig3.suptitle("{}: FSD based on PE reads".format(name_file), fontsize=14)
            ax2 = fig3.add_subplot(1, 1, 1)
            ticks = numpy.arange(1, 22)
            ticks1 = map(str, ticks)
            ticks1[len(ticks1) - 1] = ">20"
            reads = []
            reads_rel = []

            # barWidth = 0 - (len(list_to_plot) + 1) / 2 * 1. / (len(list_to_plot) + 1)
            ax2.set_xticks([], [])

            list_y = []
            label = ["duplex", "ab", "ba"]
            col = ["#FF0000", "#5FB404", "#FFBF00"]
            for i in range(len(list1)):
                x = list(numpy.arange(1, 22).astype(float))
                unique, c = numpy.unique(list1[i], return_counts=True)
                y = unique * c
                if sum(list1_o[i] > 20) > 0:
                    y[len(y) - 1] = sum(list1_o[i][list1_o[i] > 20])
                y = [y[x[idx] == unique][0] if x[idx] in unique else 0 for idx in range(len(x))]
                reads.append(y)
                reads_rel.append(list(numpy.float_(y)) / sum(numpy.concatenate(list1_o)))

                if rel_freq:
                    y = list(numpy.float_(y)) / sum(numpy.concatenate(list1_o))
                    ax2.set_ylim(0, 1.07)
                else:
                    y = y

                list_y.append(y)
                if i == 0:
                    ax2.bar(x, y, align="center", width=0.8, edgecolor="black", label=label[0], linewidth=1, alpha=1, color=col[0])
                elif i == 1:
                    ax2.bar(x, y, bottom=list_y[i - 1], align="center", width=0.8, edgecolor="black", label=label[1], linewidth=1, alpha=1, color=col[1])
                elif i == 2:
                    bars = numpy.add(list_y[0], list_y[1]).tolist()
                    ax2.bar(x, y, bottom=bars, align="center", width=0.8, edgecolor="black", label=label[2], linewidth=1, alpha=1, color=col[2])

            ax2.legend(loc='upper right', fontsize=14, frameon=True, bbox_to_anchor=(0.9, 1))

            singl = len(data_o[data_o == 1])
            last = len(data_o[data_o > 20])  # large families

            ax2.set_xticks(numpy.array(ticks))
            ax2.set_xticklabels(ticks1)
            ax2.set_xlabel("Family size", fontsize=14)
            ax2.set_ylabel(ylab, fontsize=14)
            if log_axis:
                ax2.set_yscale('log')
            ax2.grid(b=True, which="major", color="#424242", linestyle=":")
            ax2.margins(0.01, None)

            # extra information beneath the plot
            legend = "SSCS ab= \nSSCS ba= \nDCS (total)= \ntotal nr. of tags="
            plt.text(0.1, 0.09, legend, size=10, transform=plt.gcf().transFigure)

            legend = "nr. of tags\n\n{:,}\n{:,}\n{:,} ({:,})\n{:,} ({:,})".format(len(dataAB), len(dataBA),
                                                                                  len(duplTags), len(duplTags_double), (len(dataAB) + len(dataBA) + len(duplTags)),
                                                                                  (len(ab) + len(ba)))
            plt.text(0.23, 0.09, legend, size=10, transform=plt.gcf().transFigure)

            legend5 = "PE reads\n\n{:,}\n{:,}\n{:,} ({:,})\n{:,} ({:,})".format(sum(dataAB_o), sum(dataBA_o),
                                                                                sum(duplTags_o), sum(duplTags_double_o),
                                                                                (sum(dataAB_o) + sum(dataBA_o) + sum(duplTags_o)),
                                                                                (sum(ab_o) + sum(ba_o)))
            plt.text(0.38, 0.09, legend5, size=10, transform=plt.gcf().transFigure)

            legend = "rel. freq. of tags\nunique\n{:.3f}\n{:.3f}\n{:.3f}\n{:,}".format(
                float(len(dataAB)) / (len(dataAB) + len(dataBA) + len(duplTags)),
                float(len(dataBA)) / (len(dataAB) + len(dataBA) + len(duplTags)),
                float(len(duplTags)) / (len(dataAB) + len(dataBA) + len(duplTags)),
                (len(dataAB) + len(dataBA) + len(duplTags)))
            plt.text(0.54, 0.09, legend, size=10, transform=plt.gcf().transFigure)

            legend = "total\n{:.3f}\n{:.3f}\n{:.3f} ({:.3f})\n{:,}".format(float(len(dataAB)) / (len(ab) + len(ba)),
                                                                           float(len(dataBA)) / (len(ab) + len(ba)),
                                                                           float(len(duplTags)) / (len(ab) + len(ba)),
                                                                           float(len(duplTags_double)) / (len(ab) + len(ba)),
                                                                           (len(ab) + len(ba)))
            plt.text(0.64, 0.09, legend, size=10, transform=plt.gcf().transFigure)

            legend1 = "\nsingletons:\nfamily size > 20:"
            plt.text(0.1, 0.03, legend1, size=10, transform=plt.gcf().transFigure)

            legend4 = "{:,}\n{:,}".format(singl, last)
            plt.text(0.23, 0.03, legend4, size=10, transform=plt.gcf().transFigure)
            legend3 = "{:.3f}\n{:.3f}".format(float(singl) / len(data), float(last) / len(data))
            plt.text(0.64, 0.03, legend3, size=10, transform=plt.gcf().transFigure)

            legend3 = "\n\n{:,}".format(sum(data_o[data_o > 20]))
            plt.text(0.38, 0.03, legend3, size=10, transform=plt.gcf().transFigure)

            legend3 = "{:.3f}\n{:.3f}".format(float(singl) / sum(data_o), float(sum(data_o[data_o > 20])) / sum(data_o))
            plt.text(0.84, 0.03, legend3, size=10, transform=plt.gcf().transFigure)

            legend = "PE reads\nunique\n{:.3f}\n{:.3f}\n{:.3f}\n{:,}".format(
                float(sum(dataAB_o)) / (sum(dataAB_o) + sum(dataBA_o) + sum(duplTags_o)),
                float(sum(dataBA_o)) / (sum(dataAB_o) + sum(dataBA_o) + sum(duplTags_o)),
                float(sum(duplTags_o)) / (sum(dataAB_o) + sum(dataBA_o) + sum(duplTags_o)),
                (sum(dataAB_o) + sum(dataBA_o) + sum(duplTags_o)))
            plt.text(0.74, 0.09, legend, size=10, transform=plt.gcf().transFigure)

            legend = "total\n{:.3f}\n{:.3f}\n{:.3f} ({:.3f})\n{:,}".format(
                float(sum(dataAB_o)) / (sum(ab_o) + sum(ba_o)),
                float(sum(dataBA_o)) / (sum(ab_o) + sum(ba_o)),
                float(sum(duplTags_o)) / (sum(ab_o) + sum(ba_o)),
                float(sum(duplTags_double_o)) / (sum(ab_o) + sum(ba_o)), (sum(ab_o) + sum(ba_o)))
            plt.text(0.84, 0.09, legend, size=10, transform=plt.gcf().transFigure)

            pdf.savefig(fig3)
            plt.close()

            # write same information to a csv file
            count = numpy.bincount(data_o)  # original counts of family sizes

            output_file.write("\nDataset:{}{}\n".format(sep, name_file))
            output_file.write("max. family size:{}{}\n".format(sep, max(data_o)))
            output_file.write("absolute frequency:{}{}\n".format(sep, count[len(count) - 1]))
            output_file.write("relative frequency:{}{:.3f}\n\n".format(sep, float(count[len(count) - 1]) / sum(count)))

            output_file.write("median family size:{}{}\n".format(sep, numpy.median(numpy.array(data_o))))
            output_file.write("mean family size:{}{}\n\n".format(sep, numpy.mean(numpy.array(data_o))))

            output_file.write(
                "{}singletons:{}{}{}family size > 20:{}{}{}{}length of dataset:\n".format(sep, sep, sep, sep, sep, sep,
                                                                                          sep, sep))
            output_file.write(
                "{}nr. of tags{}rel. freq of tags{}rel.freq of PE reads{}nr. of tags{}rel. freq of tags{}nr. of PE reads{}rel. freq of PE reads{}total nr. of tags{}total nr. of PE reads\n".format(
                    sep, sep, sep, sep, sep, sep, sep, sep, sep))
            output_file.write("{}{}{}{}{:.3f}{}{:.3f}{}{}{}{:.3f}{}{}{}{:.3f}{}{}{}{}\n\n".format(
                name_file, sep, singl, sep, float(singl) / len(data), sep, float(singl) / sum(data_o), sep,
                last, sep, float(last) / len(data), sep, sum(data_o[data_o > 20]), sep, float(sum(data_o[data_o > 20])) / sum(data_o), sep, len(data),
                sep, sum(data_o)))

            # information for FS >= 1
            output_file.write(
                "The unique frequencies were calculated from the dataset where the tags occured only once (=ab without DCS, ba without DCS)\n"
                "Whereas the total frequencies were calculated from the whole dataset (=including the DCS).\n\n")
            output_file.write(
                "FS >= 1{}nr. of tags{}nr. of PE reads{}rel. freq of tags{}{}rel. freq of PE reads:\n".format(sep, sep,
                                                                                                              sep, sep,
                                                                                                              sep))
            output_file.write("{}{}{}unique:{}total{}unique{}total:\n".format(sep, sep, sep, sep, sep, sep))
            output_file.write("SSCS ab{}{}{}{}{}{:.3f}{}{:.3f}{}{:.3f}{}{:.3f}\n".format(
                sep, len(dataAB), sep, sum(dataAB_o), sep,
                float(len(dataAB)) / (len(dataAB) + len(dataBA) + len(duplTags)),
                sep, float(len(dataAB)) / (len(ab) + len(ba)), sep, float(sum(dataAB_o)) / (sum(dataAB_o) + sum(dataBA_o) + sum(duplTags_o)),
                sep, float(sum(dataAB_o)) / (sum(ab_o) + sum(ba_o))))
            output_file.write("SSCS ba{}{}{}{}{}{:.3f}{}{:.3f}{}{:.3f}{}{:.3f}\n".format(
                sep, len(dataBA), sep, sum(dataBA_o), sep,
                float(len(dataBA)) / (len(dataAB) + len(dataBA) + len(duplTags)),
                sep, float(len(dataBA)) / (len(ab) + len(ba)), sep,
                float(sum(dataBA_o)) / (sum(dataAB_o) + sum(dataBA_o) + sum(duplTags_o)),
                sep, float(sum(dataBA_o)) / (sum(ab_o) + sum(ba_o))))
            output_file.write(
                "DCS (total){}{} ({}){}{} ({}){}{:.3f}{}{:.3f} ({:.3f}){}{:.3f}{}{:.3f} ({:.3f})\n".format(
                    sep, len(duplTags), len(duplTags_double), sep, sum(duplTags_o), sum(duplTags_double_o), sep,
                    float(len(duplTags)) / (len(dataAB) + len(dataBA) + len(duplTags)), sep,
                    float(len(duplTags)) / (len(ab) + len(ba)), float(len(duplTags_double)) / (len(ab) + len(ba)), sep,
                    float(sum(duplTags_o)) / (sum(dataAB_o) + sum(dataBA_o) + sum(duplTags_o)), sep,
                    float(sum(duplTags_o)) / (sum(ab_o) + sum(ba_o)),
                    float(sum(duplTags_double_o)) / (sum(ab_o) + sum(ba_o))))
            output_file.write("total nr. of tags{}{}{}{}{}{}{}{}{}{}{}{}\n".format(
                sep, (len(dataAB) + len(dataBA) + len(duplTags)), sep,
                (sum(dataAB_o) + sum(dataBA_o) + sum(duplTags_o)), sep,
                (len(dataAB) + len(dataBA) + len(duplTags)), sep, (len(ab) + len(ba)), sep,
                (sum(dataAB_o) + sum(dataBA_o) + sum(duplTags_o)), sep, (sum(ab_o) + sum(ba_o))))
            # information for FS >= 3
            output_file.write(
                "\nFS >= 3{}nr. of tags{}nr. of PE reads{}rel. freq of tags{}{}rel. freq of PE reads:\n".format(sep,
                                                                                                                sep,
                                                                                                                sep,
                                                                                                                sep,
                                                                                                                sep))
            output_file.write("{}{}{}unique:{}total{}unique{}total:\n".format(sep, sep, sep, sep, sep, sep))
            output_file.write("SSCS ab{}{}{}{}{}{:.3f}{}{:.3f}{}{:.3f}{}{:.3f}\n".format(
                sep, len(dataAB_FS3), sep, sum(dataAB_FS3_o), sep,
                float(len(dataAB_FS3)) / (len(dataAB_FS3) + len(dataBA_FS3) + len(duplTags_FS3)), sep,
                float(len(dataAB_FS3)) / (len(dataBA_FS3) + len(dataBA_FS3) + duplTags_double_FS3),
                sep, float(sum(dataAB_FS3_o)) / (sum(dataAB_FS3_o) + sum(dataBA_FS3_o) + sum(duplTags_FS3_o)),
                sep, float(sum(dataAB_FS3_o)) / (sum(dataBA_FS3_o) + sum(dataBA_FS3_o) + duplTags_double_FS3_o)))
            output_file.write("SSCS ba{}{}{}{}{}{:.3f}{}{:.3f}{}{:.3f}{}{:.3f}\n".format(
                sep, len(dataBA_FS3), sep, sum(dataBA_FS3_o), sep,
                float(len(dataBA_FS3)) / (len(dataBA_FS3) + len(dataBA_FS3) + len(duplTags_FS3)),
                sep, float(len(dataBA_FS3)) / (len(dataBA_FS3) + len(dataBA_FS3) + duplTags_double_FS3),
                sep, float(sum(dataBA_FS3_o)) / (sum(dataBA_FS3_o) + sum(dataBA_FS3_o) + sum(duplTags_FS3_o)),
                sep, float(sum(dataBA_FS3_o)) / (sum(dataBA_FS3_o) + sum(dataBA_FS3_o) + duplTags_double_FS3_o)))
            output_file.write(
                "DCS (total){}{} ({}){}{} ({}){}{:.3f}{}{:.3f} ({:.3f}){}{:.3f}{}{:.3f} ({:.3f})\n".format(
                    sep, len(duplTags_FS3), duplTags_double_FS3, sep, sum(duplTags_FS3_o), duplTags_double_FS3_o, sep,
                    float(len(duplTags_FS3)) / (len(dataAB_FS3) + len(dataBA_FS3) + len(duplTags_FS3)), sep,
                    float(len(duplTags_FS3)) / (len(dataAB_FS3) + len(dataBA_FS3) + duplTags_double_FS3),
                    float(duplTags_double_FS3) / (len(dataAB_FS3) + len(dataBA_FS3) + duplTags_double_FS3),
                    sep, float(sum(duplTags_FS3_o)) / (sum(dataAB_FS3_o) + sum(dataBA_FS3_o) + sum(duplTags_FS3_o)),
                    sep,
                    float(sum(duplTags_FS3_o)) / (sum(dataAB_FS3_o) + sum(dataBA_FS3_o) + duplTags_double_FS3_o),
                    float(duplTags_double_FS3_o) / (sum(dataAB_FS3_o) + sum(dataBA_FS3_o) + duplTags_double_FS3_o)))
            output_file.write("total nr. of tags{}{}{}{}{}{}{}{}{}{}{}{}\n".format(
                sep, (len(dataAB_FS3) + len(dataBA_FS3) + len(duplTags_FS3)), sep,
                (sum(dataAB_FS3_o) + sum(dataBA_FS3_o) + sum(duplTags_FS3_o)),
                sep, (len(dataAB_FS3) + len(dataBA_FS3) + len(duplTags_FS3)), sep,
                (len(dataAB_FS3) + len(dataBA_FS3) + duplTags_double_FS3),
                sep, (sum(dataAB_FS3_o) + sum(dataBA_FS3_o) + sum(duplTags_FS3_o)), sep,
                (sum(dataAB_FS3_o) + sum(dataBA_FS3_o) + duplTags_double_FS3_o)))

            counts = [numpy.bincount(dk, minlength=22)[1:] for dk in list1]  # original counts of family sizes
            output_file.write("\nValues from family size distribution based on families\n")
            output_file.write("{}duplex{}ab{}ba{}sum\n".format(sep, sep, sep, sep))

            j = 0
            for fs in bins:
                if fs == 21:
                    fs = ">20"
                else:
                    fs = "={}".format(fs)
                output_file.write("FS{}{}".format(fs, sep))
                for n in range(3):
                    output_file.write("{}{}".format(int(counts[n][j]), sep))
                output_file.write("{}\n".format(counts[0][j] + counts[1][j] + counts[2][j]))
                j += 1
            output_file.write("sum{}".format(sep))
            for i in counts:
                output_file.write("{}{}".format(int(sum(i)), sep))
            output_file.write("{}\n".format(sum(counts[0] + counts[1] + counts[2])))

            output_file.write("\nValues from family size distribution based on PE reads\n")
            output_file.write("{}duplex{}ab{}ba{}sum\n".format(sep, sep, sep, sep))
            j = 0
            for fs in bins:
                if fs == 21:
                    fs = ">20"
                else:
                    fs = "={}".format(fs)
                output_file.write("FS{}{}".format(fs, sep))
                for n in range(3):
                    output_file.write("{}{}".format(int(reads[n][j]), sep))
                output_file.write("{}\n".format(reads[0][j] + reads[1][j] + reads[2][j]))
                j += 1
            output_file.write("sum{}".format(sep))
            for i in reads:
                output_file.write("{}{}".format(int(sum(i)), sep))
            output_file.write("{}\n".format(sum(reads[0] + reads[1] + reads[2])))

    print("Files successfully created!")


if __name__ == '__main__':
    sys.exit(compare_read_families(sys.argv))
