#!/usr/bin/env python

# Family size distribution of SSCSs
#
# Author: Monika Heinzl, Johannes-Kepler University Linz (Austria)
# Contact: monika.heinzl@edumail.at
#
# Takes at least one TABULAR file with tags before the alignment to the SSCS, but up to 4 files can be provided, as input.
# The program produces a plot which shows the distribution of family sizes of the all SSCSs from the input files and
# a CSV file with the data of the plot, as well as a TXT file with all tags of the DCS and their family sizes.
# If only one file is provided, then a family size distribution, which is separated after SSCSs without a partner and DCSs, is produced.
# Whereas a family size distribution with multiple data in one plot is produced, when more than one file (up to 4) is given.

# USAGE: python FSD_Galaxy_1.4_commandLine_FINAL.py --inputFile1 filename --inputName1 filename --inputFile2 filename2 --inputName2 filename2 --inputFile3 filename3 --inputName3 filename3 --inputFile4 filename4 --inputName4 filename4 --sep "characterWhichSeparatesCSVFile" --output_csv outptufile_name_csv --output_pdf outptufile_name_pdf

import argparse
import sys

import matplotlib.pyplot as plt
import numpy
from matplotlib.backends.backend_pdf import PdfPages


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
    parser.add_argument('--sep', default=",", help='Separator in the csv file.')
    parser.add_argument('--output_pdf', default="data.pdf", type=str, help='Name of the pdf file.')
    parser.add_argument('--output_csv', default="data.csv", type=str, help='Name of the csv file.')
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

    title_file = args.output_csv
    title_file2 = args.output_pdf
    sep = args.sep

    if type(sep) is not str or len(sep) > 1:
        print("Error: --sep must be a single character.")
        exit(1)

    plt.rc('figure', figsize=(11.69, 8.27))  # A4 format
    plt.rcParams['patch.edgecolor'] = "black"
    plt.rcParams['axes.facecolor'] = "E0E0E0"  # grey background color
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    list_to_plot = []
    label = []
    data_array_list = []
    with open(title_file, "w") as output_file, PdfPages(title_file2) as pdf:
        fig = plt.figure()
        plt.subplots_adjust(bottom=0.25)
        if firstFile != str(None):
            file1 = readFileReferenceFree(firstFile)
            integers = numpy.array(file1[:, 0]).astype(int)  # keep original family sizes

            # for plot: replace all big family sizes by 22
            data1 = numpy.array(file1[:, 0]).astype(int)
            bigFamilies = numpy.where(data1 > 20)[0]
            data1[bigFamilies] = 22

            name1 = name1.split(".tabular")[0]
            list_to_plot.append(data1)
            label.append(name1)
            data_array_list.append(file1)

            legend = "\n\n\n{}".format(name1)
            plt.text(0.1, 0.11, legend, size=12, transform=plt.gcf().transFigure)
            legend1 = "singletons:\nabsolute nr.\n{:,}".format(numpy.bincount(data1)[1])
            plt.text(0.4, 0.11, legend1, size=12, transform=plt.gcf().transFigure)

            legend3 = "rel. freq\n{:.3f}".format(float(numpy.bincount(data1)[1]) / len(data1))
            plt.text(0.5, 0.11, legend3, size=12, transform=plt.gcf().transFigure)

            legend4 = "family size > 20:\nabsolute nr.\n{:,}".format(numpy.bincount(data1)[len(numpy.bincount(data1)) - 1].astype(int))
            plt.text(0.6, 0.11, legend4, size=12, transform=plt.gcf().transFigure)

            legend5 = "rel. freq\n{:.3f}".format(float(numpy.bincount(data1)[len(numpy.bincount(data1)) - 1]) / len(data1))
            plt.text(0.7, 0.11, legend5, size=12, transform=plt.gcf().transFigure)

            legend6 = "total length\n{:,}".format(len(data1))
            plt.text(0.8, 0.11, legend6, size=12, transform=plt.gcf().transFigure)

        if secondFile != str(None):
            file2 = readFileReferenceFree(secondFile)
            data2 = numpy.asarray(file2[:, 0]).astype(int)
            bigFamilies2 = numpy.where(data2 > 20)[0]
            data2[bigFamilies2] = 22

            list_to_plot.append(data2)
            name2 = name2.split(".tabular")[0]
            label.append(name2)
            data_array_list.append(file2)

            plt.text(0.1, 0.09, name2, size=12, transform=plt.gcf().transFigure)

            legend1 = "{:,}".format(numpy.bincount(data2)[1])
            plt.text(0.4, 0.09, legend1, size=12, transform=plt.gcf().transFigure)

            legend3 = "{:.3f}".format(float(numpy.bincount(data2)[1]) / len(data2))
            plt.text(0.5, 0.09, legend3, size=12, transform=plt.gcf().transFigure)

            legend4 = "{:,}".format(numpy.bincount(data2)[len(numpy.bincount(data2)) - 1].astype(int))
            plt.text(0.6, 0.09, legend4, size=12, transform=plt.gcf().transFigure)

            legend5 = "{:.3f}".format(float(numpy.bincount(data2)[len(numpy.bincount(data2)) - 1]) / len(data2))
            plt.text(0.7, 0.09, legend5, size=12, transform=plt.gcf().transFigure)

            legend6 = "{:,}".format(len(data2))
            plt.text(0.8, 0.09, legend6, size=12, transform=plt.gcf().transFigure)

        if thirdFile != str(None):
            file3 = readFileReferenceFree(thirdFile)

            data3 = numpy.asarray(file3[:, 0]).astype(int)
            bigFamilies3 = numpy.where(data3 > 20)[0]
            data3[bigFamilies3] = 22

            list_to_plot.append(data3)
            name3 = name3.split(".tabular")[0]
            label.append(name3)
            data_array_list.append(file3)

            plt.text(0.1, 0.07, name3, size=12, transform=plt.gcf().transFigure)

            legend1 = "{:,}".format(numpy.bincount(data3)[1])
            plt.text(0.4, 0.07, legend1, size=12, transform=plt.gcf().transFigure)

            legend3 = "{:.3f}".format(float(numpy.bincount(data3)[1]) / len(data3))
            plt.text(0.5, 0.07, legend3, size=12, transform=plt.gcf().transFigure)

            legend4 = "{:,}".format(numpy.bincount(data3)[len(numpy.bincount(data3)) - 1].astype(int))
            plt.text(0.6, 0.07, legend4, size=12, transform=plt.gcf().transFigure)

            legend5 = "{:.3f}".format(float(numpy.bincount(data3)[len(numpy.bincount(data3)) - 1]) / len(data3))
            plt.text(0.7, 0.07, legend5, size=12, transform=plt.gcf().transFigure)

            legend6 = "{:,}".format(len(data3))
            plt.text(0.8, 0.07, legend6, size=12, transform=plt.gcf().transFigure)

        if fourthFile != str(None):
            file4 = readFileReferenceFree(fourthFile)

            data4 = numpy.asarray(file4[:, 0]).astype(int)

            bigFamilies4 = numpy.where(data4 > 20)[0]
            data4[bigFamilies4] = 22

            list_to_plot.append(data4)
            name4 = name4.split(".tabular")[0]
            label.append(name4)
            data_array_list.append(file4)

            plt.text(0.1, 0.05, name4, size=12, transform=plt.gcf().transFigure)

            legend1 = "{:,}".format(numpy.bincount(data4)[1])
            plt.text(0.4, 0.05, legend1, size=12, transform=plt.gcf().transFigure)

            legend4 = "{:.3f}".format(float(numpy.bincount(data4)[1]) / len(data4))
            plt.text(0.5, 0.05, legend4, size=12, transform=plt.gcf().transFigure)

            legend4 = "{:,}".format(numpy.bincount(data4)[len(numpy.bincount(data4)) - 1].astype(int))
            plt.text(0.6, 0.05, legend4, size=12, transform=plt.gcf().transFigure)

            legend5 = "{:.3f}".format(float(numpy.bincount(data4)[len(numpy.bincount(data4)) - 1]) / len(data4))
            plt.text(0.7, 0.05, legend5, size=12, transform=plt.gcf().transFigure)

            legend6 = "{:,}".format(len(data4))
            plt.text(0.8, 0.05, legend6, size=12, transform=plt.gcf().transFigure)

        maximumX = numpy.amax(numpy.concatenate(list_to_plot))
        minimumX = numpy.amin(numpy.concatenate(list_to_plot))

        counts = plt.hist(list_to_plot, bins=range(minimumX, maximumX + 1), stacked=False, edgecolor="black",
                          linewidth=1, label=label, align="left", alpha=0.7, rwidth=0.8)

        ticks = numpy.arange(minimumX - 1, maximumX, 1)
        ticks1 = map(str, ticks)
        ticks1[len(ticks1) - 1] = ">20"
        plt.xticks(numpy.array(ticks), ticks1)

        plt.legend(loc='upper right', fontsize=14, frameon=True, bbox_to_anchor=(0.9, 1))
        # plt.title("Family Size Distribution", fontsize=14)
        plt.xlabel("Family size", fontsize=14)
        plt.ylabel("Absolute Frequency", fontsize=14)
        plt.margins(0.01, None)
        plt.grid(b=True, which="major", color="#424242", linestyle=":")
        pdf.savefig(fig)
        plt.close()

        # write data to CSV file
        output_file.write("Values from family size distribution with all datasets\n")
        output_file.write("\nFamily size")
        for i in label:
            output_file.write("{}{}".format(sep, i))
        # output_file.write("{}sum".format(sep))
        output_file.write("\n")
        j = 0
        for fs in counts[1][0:len(counts[1]) - 1]:
            if fs == 21:
                fs = ">20"
            else:
                fs = "={}".format(fs)
            output_file.write("FS{}{}".format(fs, sep))
            if len(label) == 1:
                output_file.write("{}{}".format(int(counts[0][j]), sep))
            else:
                for n in range(len(label)):
                    output_file.write("{}{}".format(int(counts[0][n][j]), sep))
            output_file.write("\n")
            j += 1
        output_file.write("sum{}".format(sep))
        values_for_sum = []
        if len(label) == 1:
            output_file.write("{}{}".format(int(sum(counts[0])), sep))
            values_for_sum.append(int(sum(counts[0])))
        else:
            for i in counts[0]:
                output_file.write("{}{}".format(int(sum(i)), sep))
                values_for_sum.append(int(sum(i)))

        output_file.write("{}\n".format(sum(values_for_sum)))

        # Family size distribution after DCS and SSCS
        for dataset, data, name_file in zip(list_to_plot, data_array_list, label):
            maximumX = numpy.amax(dataset)
            minimumX = numpy.amin(dataset)

            tags = numpy.array(data[:, 2])
            seq = numpy.array(data[:, 1])
            data = numpy.array(dataset)

            # find all unique tags and get the indices for ALL tags, but only once
            u, index_unique, c = numpy.unique(numpy.array(seq), return_counts=True, return_index=True)
            d = u[c > 1]

            # get family sizes, tag for duplicates
            duplTags_double = data[numpy.in1d(seq, d)]
            duplTags = duplTags_double[0::2]  # ab of DCS
            duplTagsBA = duplTags_double[1::2]  # ba of DCS

            # duplTags_double_tag = tags[numpy.in1d(seq, d)]
            # duplTags_double_seq = seq[numpy.in1d(seq, d)]

            # get family sizes for SSCS with no partner
            ab = numpy.where(tags == "ab")[0]
            abSeq = seq[ab]
            ab = data[ab]
            ba = numpy.where(tags == "ba")[0]
            baSeq = seq[ba]
            ba = data[ba]

            dataAB = ab[numpy.in1d(abSeq, d, invert=True)]
            dataBA = ba[numpy.in1d(baSeq, d, invert=True)]

            list1 = [duplTags_double, dataAB, dataBA]  # list for plotting

            # information for family size >= 3
            dataAB_FS3 = dataAB[dataAB >= 3]
            dataBA_FS3 = dataBA[dataBA >= 3]
            ab_FS3 = ab[ab >= 3]
            ba_FS3 = ba[ba >= 3]

            duplTags_FS3 = duplTags[(duplTags >= 3) & (duplTagsBA >= 3)]  # ab+ba with FS>=3
            duplTags_FS3_BA = duplTagsBA[(duplTags >= 3) & (duplTagsBA >= 3)]  # ba+ab with FS>=3
            duplTags_double_FS3 = len(duplTags_FS3) + len(duplTags_FS3_BA)  # both ab and ba strands with FS>=3

            fig = plt.figure()

            plt.subplots_adjust(bottom=0.3)
            counts = plt.hist(list1, bins=range(minimumX, maximumX + 1), stacked=True, label=["duplex", "ab", "ba"], edgecolor="black", linewidth=1, align="left", color=["#FF0000", "#5FB404", "#FFBF00"])
            # tick labels of x axis
            ticks = numpy.arange(minimumX - 1, maximumX, 1)
            ticks1 = map(str, ticks)
            ticks1[len(ticks1) - 1] = ">20"
            plt.xticks(numpy.array(ticks), ticks1)
            singl = counts[0][2][0]  # singletons
            last = counts[0][2][len(counts[0][0]) - 1]  # large families

            plt.legend(loc='upper right', fontsize=14, bbox_to_anchor=(0.9, 1), frameon=True)
            # plt.title(name1, fontsize=14)
            plt.xlabel("Family size", fontsize=14)
            plt.ylabel("Absolute Frequency", fontsize=14)
            plt.margins(0.01, None)
            plt.grid(b=True, which="major", color="#424242", linestyle=":")

            # extra information beneath the plot
            legend = "SSCS ab= \nSSCS ba= \nDCS (total)= \nlength of dataset="
            plt.text(0.1, 0.09, legend, size=12, transform=plt.gcf().transFigure)

            legend = "absolute numbers\n\n{:,}\n{:,}\n{:,} ({:,})\n{:,}".format(len(dataAB), len(dataBA), len(duplTags), len(duplTags_double), (len(dataAB) + len(dataBA) + len(duplTags)))
            plt.text(0.35, 0.09, legend, size=12, transform=plt.gcf().transFigure)

            legend = "relative frequencies\nunique\n{:.3f}\n{:.3f}\n{:.3f}\n{:,}".format(float(len(dataAB)) / (len(dataAB) + len(dataBA) + len(duplTags)), float(len(dataBA)) / (len(dataAB) + len(dataBA) + len(duplTags)), float(len(duplTags)) / (len(dataAB) + len(dataBA) + len(duplTags)), (len(dataAB) + len(dataBA) + len(duplTags)))
            plt.text(0.54, 0.09, legend, size=12, transform=plt.gcf().transFigure)

            legend = "total\n{:.3f}\n{:.3f}\n{:.3f} ({:.3f})\n{:,}".format(float(len(dataAB)) / (len(ab) + len(ba)), float(len(dataBA)) / (len(ab) + len(ba)), float(len(duplTags)) / (len(ab) + len(ba)), float(len(duplTags_double)) / (len(ab) + len(ba)), (len(ab) + len(ba)))
            plt.text(0.64, 0.09, legend, size=12, transform=plt.gcf().transFigure)

            legend1 = "\nsingletons:\nfamily size > 20:"
            plt.text(0.1, 0.03, legend1, size=12, transform=plt.gcf().transFigure)

            legend4 = "{:,}\n{:,}".format(singl.astype(int), last.astype(int))
            plt.text(0.35, 0.03, legend4, size=12, transform=plt.gcf().transFigure)

            legend3 = "{:.3f}\n{:.3f}".format(singl / len(data), last / len(data))
            plt.text(0.54, 0.03, legend3, size=12, transform=plt.gcf().transFigure)

            pdf.savefig(fig)
            plt.close()

            # write same information to a csv file
            count = numpy.bincount(integers)  # original counts of family sizes
            output_file.write("\nDataset:{}{}\n".format(sep, name_file))
            output_file.write("max. family size:{}{}\n".format(sep, max(integers)))
            output_file.write("absolute frequency:{}{}\n".format(sep, count[len(count) - 1]))
            output_file.write("relative frequency:{}{:.3f}\n\n".format(sep, float(count[len(count) - 1]) / sum(count)))

            output_file.write("{}singletons:{}{}family size > 20:\n".format(sep, sep, sep))
            output_file.write("{}absolute nr.{}rel. freq{}absolute nr.{}rel. freq{}total length\n".format(sep, sep, sep, sep, sep))
            output_file.write("{}{}{}{}{:.3f}{}{}{}{:.3f}{}{}\n\n".format(name_file, sep, singl.astype(int), sep, singl / len(data), sep, last.astype(int), sep, last / len(data), sep, len(data)))

            # information for FS >= 1
            output_file.write("The unique frequencies were calculated from the dataset where the tags occured only once (=ab without DCS, ba without DCS)\nWhereas the total frequencies were calculated from the whole dataset (=including the DCS).\n\n")
            output_file.write("FS >= 1{}{}unique:{}total:\n".format(sep, sep, sep))
            output_file.write("nr./rel. freq of ab={}{}{}{:.3f}{}{:.3f}\n".format(sep, len(dataAB), sep, float(len(dataAB)) / (len(dataAB) + len(dataBA) + len( duplTags)), sep, float(len(dataAB)) / (len(ab) + len(ba))))
            output_file.write("nr./rel. freq of ba={}{}{}{:.3f}{}{:.3f}\n".format(sep, len(dataBA), sep, float(len(dataBA)) / (len(dataBA) + len(dataBA) + len(duplTags)), sep, float(len(dataBA)) / (len(ba) + len(ba))))
            output_file.write("nr./rel. freq of DCS (total)={}{} ({}){}{:.3f}{}{:.3f} ({:.3f})\n".format(sep, len(duplTags), len(duplTags_double), sep, float(len(duplTags)) / (len(dataAB) + len(dataBA) + len(duplTags)), sep, float(len(duplTags)) / ( len(ab) + len(ba)), float(len(duplTags_double)) / (len(ab) + len(ba))))
            output_file.write("length of dataset={}{}{}{}{}{}\n".format(sep, (len(dataAB) + len(dataBA) + len(duplTags)), sep, (len(dataAB) + len(dataBA) + len(duplTags)), sep, (len(ab) + len(ba))))
            # information for FS >= 3
            output_file.write("FS >= 3{}{}unique:{}total:\n".format(sep, sep, sep))
            output_file.write("nr./rel. freq of ab={}{}{}{:.3f}{}{:.3f}\n".format(sep, len(dataAB_FS3), sep, float(len(dataAB_FS3)) / (len(dataAB_FS3) + len(dataBA_FS3) + len(duplTags_FS3)), sep, float(len(dataAB_FS3)) / (len(ab_FS3) + len(ba_FS3))))
            output_file.write("nr./rel. freq of ba={}{}{}{:.3f}{}{:.3f}\n".format(sep, len(dataBA_FS3), sep, float(len(dataBA_FS3)) / (len(dataBA_FS3) + len(dataBA_FS3) + len(duplTags_FS3)), sep, float(len(dataBA_FS3)) / (len(ba_FS3) + len(ba_FS3))))
            output_file.write("nr./rel. freq of DCS (total)={}{} ({}){}{:.3f}{}{:.3f} ({:.3f})\n".format(sep, len(duplTags_FS3), duplTags_double_FS3, sep, float(len( duplTags_FS3)) / (len(dataBA_FS3) + len(duplTags_FS3)), sep, float(len(duplTags_FS3)) / (len(ab_FS3) + len(ba_FS3)), float(duplTags_double_FS3) / (len(ab_FS3) + len(ba_FS3))))
            output_file.write("length of dataset={}{}{}{}{}{}\n".format(sep, (len(dataAB_FS3) + len(dataBA_FS3) + len(duplTags_FS3)), sep, (len(dataAB_FS3) + len(dataBA_FS3) + len(duplTags_FS3)), sep, (len(ab_FS3) + len(ba_FS3))))

            output_file.write("\nValues from family size distribution\n")
            output_file.write("{}duplex{}ab{}ba{}sum\n".format(sep, sep, sep, sep))
            for dx, ab, ba, fs in zip(counts[0][0], counts[0][1], counts[0][2], counts[1]):
                if fs == 21:
                    fs = ">20"
                else:
                    fs = "={}".format(fs)
                ab1 = ab - dx
                ba1 = ba - ab
                output_file.write("FS{}{}{}{}{}{}{}{}{}\n".format(fs, sep, int(dx), sep, int(ab1), sep, int(ba1), sep, int(ba)))

    print("Files successfully created!")


if __name__ == '__main__':
    sys.exit(compare_read_families(sys.argv))
