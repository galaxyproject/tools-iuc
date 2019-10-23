#!/usr/bin/env python
# Family size distribution of DCS from various steps of the Galaxy pipeline
#
# Author: Monika Heinzl & Gundula Povysil, Johannes-Kepler University Linz (Austria)
# Contact: monika.heinzl@edumail.at
#
# Takes a TXT file with tags of reads that were aligned to certain regions of the reference genome (optional),
# a TABULAR file with tags before the alignment to the SSCS, a FASTA file with reads that were part of the DCS and
# a FASTA file with tags after trimming as input (optional).
# The program produces a plot which shows the distribution of family sizes of the DCS from the input files and
# a CSV file with the data of the plot.
# USAGE: python FSD before vs after_no_refF1.3_FINAL.py --inputFile_SSCS filenameSSCS --inputName1 filenameSSCS --makeDCS filenameMakeDCS --afterTrimming filenameAfterTrimming --alignedTags DCSbamFile
# --output_tabular outputfile_name_tabular --output_pdf outputfile_name_pdf

import argparse
import re
import sys
from collections import Counter

import matplotlib.pyplot as plt
import numpy
import pysam
from Bio import SeqIO
from matplotlib.backends.backend_pdf import PdfPages

plt.switch_backend('agg')


def readFileReferenceFree(file, delim):
    with open(file, 'r') as dest_f:
        data_array = numpy.genfromtxt(dest_f, skip_header=0, delimiter=delim, comments='#', dtype='string')
        return data_array


def readFasta(file):
    tag_consensus = []
    fs_consensus = []
    with open(file, "r") as consFile:
        for record in SeqIO.parse(consFile, "fasta"):
            tag_consensus.append(record.id)
            line = record.description
            a, b = line.split(" ")
            fs1, fs2 = b.split("-")
            fs_consensus.extend([fs1, fs2])
    fs_consensus = numpy.array(fs_consensus).astype(int)
    return(tag_consensus, fs_consensus)


def make_argparser():
    parser = argparse.ArgumentParser(description='Analysis of read loss in duplex sequencing data')
    parser.add_argument('--inputFile_SSCS',
                        help='Tabular File with three columns: ab or ba, tag and family size.')
    parser.add_argument('--inputName1')
    parser.add_argument('--makeDCS',
                        help='FASTA File with information about tag and family size in the header.')
    parser.add_argument('--afterTrimming', default=None,
                        help='FASTA File with information about tag and family size in the header.')
    parser.add_argument('--bamFile',
                        help='BAM file with aligned reads.')
    parser.add_argument('--output_pdf', default="data.pdf", type=str,
                        help='Name of the pdf and tabular file.')
    parser.add_argument('--output_tabular', default="data.tabular", type=str,
                        help='Name of the pdf and tabular file.')
    return parser


def compare_read_families_read_loss(argv):
    parser = make_argparser()
    args = parser.parse_args(argv[1:])
    SSCS_file = args.inputFile_SSCS
    SSCS_file_name = args.inputName1
    makeConsensus = args.makeDCS
    afterTrimming = args.afterTrimming
    ref_genome = args.bamFile
    title_file = args.output_tabular
    title_file2 = args.output_pdf
    sep = "\t"

    with open(title_file, "w") as output_file, PdfPages(title_file2) as pdf:
        # PLOT
        plt.rc('figure', figsize=(11.69, 8.27))  # A4 format
        plt.rcParams['axes.facecolor'] = "E0E0E0"  # grey background color
        plt.rcParams['xtick.labelsize'] = 14
        plt.rcParams['ytick.labelsize'] = 14
        plt.rcParams['patch.edgecolor'] = "black"
        fig = plt.figure()
        plt.subplots_adjust(bottom=0.3)

        list1 = []
        colors = []
        labels = []

# data with tags of SSCS
        data_array = readFileReferenceFree(SSCS_file, "\t")
        seq = numpy.array(data_array[:, 1])
        tags = numpy.array(data_array[:, 2])
        quant = numpy.array(data_array[:, 0]).astype(int)

        # split data with all tags of SSCS after ab and ba strands
        all_ab = seq[numpy.where(tags == "ab")[0]]
        all_ba = seq[numpy.where(tags == "ba")[0]]
        quant_ab_sscs = quant[numpy.where(tags == "ab")[0]]
        quant_ba_sscs = quant[numpy.where(tags == "ba")[0]]

        seqDic_ab = dict(zip(all_ab, quant_ab_sscs))
        seqDic_ba = dict(zip(all_ba, quant_ba_sscs))

        # get tags of the SSCS which form a DCS
        # group large family sizes
        bigFamilies = numpy.where(quant > 20)[0]
        quant[bigFamilies] = 22
        maximumX = numpy.amax(quant)

        # find all unique tags and get the indices for ALL tags (ab AND ba)
        u, index_unique, c = numpy.unique(numpy.array(seq), return_counts=True, return_index=True)
        d = u[c > 1]

        # get family sizes, tag for the duplicates
        duplTags_double = quant[numpy.in1d(seq, d)]
        list1.append(duplTags_double)
        colors.append("#0000FF")
        labels.append("before SSCS building")

        duplTags = duplTags_double[0::2]  # ab of DCS
        duplTagsBA = duplTags_double[1::2]  # ba of DCS

        d2 = d[(duplTags >= 3) & (duplTagsBA >= 3)]  # ab and ba FS>=3

        # all SSCSs FS>=3
        seq_unique, seqUnique_index = numpy.unique(seq, return_index=True)
        seq_unique_FS = quant[seqUnique_index]
        seq_unique_FS3 = seq_unique_FS[seq_unique_FS >= 3]

        legend1 = "\ntotal nr. of tags (unique, FS>=1):\nDCS (before SSCS building, FS>=1):\ntotal nr. of tags (unique, FS>=3):\nDCS (before SSCS building, FS>=3):"
        legend2 = "total numbers * \n{:,}\n{:,}\n{:,}\n{:,}".format(len(seq_unique_FS), len(duplTags),
                                                                    len(seq_unique_FS3), len(d2))
        plt.text(0.55, 0.14, legend1, size=11, transform=plt.gcf().transFigure)
        plt.text(0.88, 0.14, legend2, size=11, transform=plt.gcf().transFigure)

        # data make DCS
        tag_consensus, fs_consensus = readFasta(makeConsensus)
        # group large family sizes in the plot of fasta files
        bigFamilies = numpy.where(fs_consensus > 20)[0]
        fs_consensus[bigFamilies] = 22
        list1.append(fs_consensus)
        colors.append("#298A08")
        labels.append("after DCS building")
        legend3 = "after DCS building:"
        legend4 = "{:,}".format(len(tag_consensus))
        plt.text(0.55, 0.11, legend3, size=11, transform=plt.gcf().transFigure)
        plt.text(0.88, 0.11, legend4, size=11, transform=plt.gcf().transFigure)

        # data after trimming
        if afterTrimming is not None:
            tag_trimming, fs_trimming = readFasta(afterTrimming)
            bigFamilies = numpy.where(fs_trimming > 20)[0]
            fs_trimming[bigFamilies] = 22
            list1.append(fs_trimming)
            colors.append("#DF0101")
            labels.append("after trimming")
            legend5 = "after trimming:"
            legend6 = "{:,}".format(len(tag_trimming))
            plt.text(0.55, 0.09, legend5, size=11, transform=plt.gcf().transFigure)
            plt.text(0.88, 0.09, legend6, size=11, transform=plt.gcf().transFigure)

# data of tags aligned to reference genome
        if ref_genome is not None:
            pysam.index(ref_genome)
            bam = pysam.AlignmentFile(ref_genome, "rb")
            seq_mut = []
            for read in bam.fetch():
                if not read.is_unmapped:
                    if '_' in read.query_name:
                        tags = read.query_name.split('_')[0]
                    else:
                        tags = read.query_name
                    seq_mut.append(tags)

            # use only unique tags that were alignment to the reference genome
            seq_mut = numpy.array(seq_mut)
            seq_mut, seqMut_index = numpy.unique(seq_mut, return_index=True)
            # get family sizes for each tag in the BAM file
            quant_ab = []
            quant_ba = []
            for i in seq_mut:
                quant_ab.append(seqDic_ab.get(i))
                quant_ba.append(seqDic_ba.get(i))

            quant_ab_ref = numpy.array(quant_ab)
            quant_ba_ref = numpy.array(quant_ba)
            quant_all_ref = numpy.concatenate((quant_ab_ref, quant_ba_ref))
            bigFamilies = numpy.where(quant_all_ref > 20)[0]  # group large family sizes
            quant_all_ref[bigFamilies] = 22

            list1.append(quant_all_ref)
            colors.append("#04cec7")
            labels.append("after alignment\nto reference")
            legend7 = "after alignment to reference:"
            length_DCS_ref = len(quant_ba_ref)  # count of duplex tags that were aligned to reference genome
            legend8 = "{:,}".format(length_DCS_ref)
            plt.text(0.55, 0.07, legend7, size=11, transform=plt.gcf().transFigure)
            plt.text(0.88, 0.07, legend8, size=11, transform=plt.gcf().transFigure)

        counts = plt.hist(list1, bins=range(-1, maximumX + 1), stacked=False, label=labels, color=colors,
                          align="left", alpha=1, edgecolor="black", linewidth=1)
        ticks = numpy.arange(0, maximumX, 1)
        ticks1 = map(str, ticks)
        ticks1[len(ticks1) - 1] = ">20"
        plt.xticks(numpy.array(ticks), ticks1)
        if ref_genome is not None:
            count = numpy.array([v for k, v in sorted(Counter(quant_ab_ref).iteritems())])  # count all family sizes from all ab strands

            legend = "max. family size:\nabsolute frequency:\nrelative frequency:\n\ntotal nr. of reads:\n(before SSCS building)"
            plt.text(0.1, 0.085, legend, size=11, transform=plt.gcf().transFigure)

            legend = "AB\n{}\n{}\n{:.5f}\n\n{:,}" \
                .format(max(quant_ab_ref), count[len(count) - 1], float(count[len(count) - 1]) / sum(count),
                        sum(numpy.array(data_array[:, 0]).astype(int)))
            plt.text(0.35, 0.105, legend, size=11, transform=plt.gcf().transFigure)

            count2 = numpy.array(
                [v for k, v in sorted(Counter(quant_ba_ref).iteritems())])  # count all family sizes from all ba strands
            legend = "BA\n{}\n{}\n{:.5f}" \
                .format(max(quant_ba_ref), count2[len(count2) - 1], float(count2[len(count2) - 1]) / sum(count2))
            plt.text(0.45, 0.1475, legend, size=11, transform=plt.gcf().transFigure)

        legend4 = "* In the plot, the family sizes of ab and ba strands and of both duplex tags were used.\nWhereas the total numbers indicate only the single count of the formed duplex tags."
        plt.text(0.1, 0.02, legend4, size=11, transform=plt.gcf().transFigure)

        plt.legend(loc='upper right', fontsize=14, bbox_to_anchor=(0.9, 1), frameon=True)
        plt.title("Family size distribution of tags from various steps of the Du Novo pipeline", fontsize=14)
        plt.xlabel("Family size", fontsize=14)
        plt.ylabel("Absolute Frequency", fontsize=14)
        plt.grid(b=True, which="major", color="#424242", linestyle=":")
        plt.margins(0.01, None)

        pdf.savefig(fig, bbox_inch="tight")
        plt.close()

    # write information about plot into a csv file
        output_file.write("Dataset:{}{}\n".format(sep, SSCS_file_name))
        if ref_genome != str(None):
            output_file.write("{}AB{}BA\n".format(sep, sep))
            output_file.write("max. family size:{}{}{}{}\n".format(sep, max(quant_ab_ref), sep, max(quant_ba_ref)))
            output_file.write(
                "absolute frequency:{}{}{}{}\n".format(sep, count[len(count) - 1], sep, count2[len(count2) - 1]))
            output_file.write(
                "relative frequency:{}{:.3f}{}{:.3f}\n\n".format(sep, float(count[len(count) - 1]) / sum(count), sep,
                                                                 float(count2[len(count2) - 1]) / sum(count2)))

        output_file.write("\ntotal nr. of reads before SSCS building{}{}\n".format(sep, sum(numpy.array(data_array[:, 0]).astype(int))))
        output_file.write("\n\nValues from family size distribution\n")

        if afterTrimming == str(None) and ref_genome == str(None):
            if afterTrimming == str(None):
                output_file.write("{}before SSCS building{}after DCS building\n".format(sep, sep))
            elif ref_genome == str(None):
                output_file.write("{}before SSCS building{}atfer DCS building\n".format(sep, sep))

            for fs, sscs, dcs in zip(counts[1][2:len(counts[1])], counts[0][0][2:len(counts[0][0])], counts[0][1][2:len(counts[0][1])]):
                if fs == 21:
                    fs = ">20"
                else:
                    fs = "={}".format(fs)
                output_file.write("FS{}{}{}{}{}\n".format(fs, sep, int(sscs), sep, int(dcs)))
            output_file.write("sum{}{}{}{}\n".format(sep, int(sum(counts[0][0])), sep, int(sum(counts[0][1]))))

        elif afterTrimming == str(None) or ref_genome == str(None):
            if afterTrimming == str(None):
                output_file.write("{}before SSCS building{}after DCS building{}after alignment to reference\n".format(sep, sep, sep))
            elif ref_genome == str(None):
                output_file.write("{}before SSCS building{}atfer DCS building{}after trimming\n".format(sep, sep, sep))

            for fs, sscs, dcs, reference in zip(counts[1][2:len(counts[1])], counts[0][0][2:len(counts[0][0])], counts[0][1][2:len(counts[0][1])], counts[0][2][2:len(counts[0][2])]):
                if fs == 21:
                    fs = ">20"
                else:
                    fs = "={}".format(fs)
                output_file.write("FS{}{}{}{}{}{}{}\n".format(fs, sep, int(sscs), sep, int(dcs), sep, int(reference)))
            output_file.write("sum{}{}{}{}{}{}\n".format(sep, int(sum(counts[0][0])), sep, int(sum(counts[0][1])), sep, int(sum(counts[0][2]))))

        else:
            output_file.write("{}before SSCS building{}after DCS building{}after trimming{}after alignment to reference\n".format(sep, sep, sep, sep))
            for fs, sscs, dcs, trim, reference in zip(counts[1][2:len(counts[1])], counts[0][0][2:len(counts[0][0])], counts[0][1][2:len(counts[0][1])], counts[0][2][2:len(counts[0][2])], counts[0][3][2:len(counts[0][3])]):
                if fs == 21:
                    fs = ">20"
                else:
                    fs = "={}".format(fs)
                output_file.write("FS{}{}{}{}{}{}{}{}{}\n".format(fs, sep, int(sscs), sep, int(dcs), sep, int(trim), sep, int(reference)))
            output_file.write("sum{}{}{}{}{}{}{}{}\n".format(sep, int(sum(counts[0][0])), sep, int(sum(counts[0][1])), sep, int(sum(counts[0][2])), sep, int(sum(counts[0][3]))))

        output_file.write("\n\nIn the plot, the family sizes of ab and ba strands and of both duplex tags were used.\nWhereas the total numbers indicate only the single count of the formed duplex tags.\n")
        output_file.write("total nr. of tags (unique, FS>=1){}{}\n".format(sep, len(seq_unique_FS)))
        output_file.write("DCS (before SSCS building, FS>=1){}{}\n".format(sep, len(duplTags)))
        output_file.write("total nr. of tags (unique, FS>=3){}{}\n".format(sep, len(seq_unique_FS3)))
        output_file.write("DCS (before SSCS building, FS>=3){}{}\n".format(sep, len(d2)))
        output_file.write("after DCS building{}{}\n".format(sep, len(tag_consensus)))
        if afterTrimming != str(None):
            output_file.write("after trimming{}{}\n".format(sep, len(tag_trimming)))
        if ref_genome is not None:
            output_file.write("after alignment to reference{}{}\n".format(sep, length_DCS_ref))

        print("Files successfully created!")


if __name__ == '__main__':
    sys.exit(compare_read_families_read_loss(sys.argv))
