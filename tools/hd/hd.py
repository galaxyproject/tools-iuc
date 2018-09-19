#!/usr/bin/env python

# Hamming distance analysis of SSCSs
#
# Author: Monika Heinzl, Johannes-Kepler University Linz (Austria)
# Contact: monika.heinzl@edumail.at
#
# Takes at least one TABULAR file with tags before the alignment to the SSCS and optionally a second TABULAR file as input.
# The program produces a plot which shows a histogram of Hamming distances separated after family sizes,
# a family size distribution separated after Hamming distances for all (sample_size=0) or a given sample of SSCSs or SSCSs, which form a DCS.
# In additon, the tool produces HD and FSD plots for the difference between the HDs of both parts of the tags and for the chimeric reads
# and finally a CSV file with the data of the plots.
# It is also possible to perform the HD analysis with shortened tags with given sizes as input.
# The tool can run on a certain number of processors, which can be defined by the user.

# USAGE: python hd.py --inputFile filename --inputName1 filename --inputFile2 filename2 --inputName2 filename2 --sample_size int/0 --sep "characterWhichSeparatesCSVFile" /
#        --only_DCS True --FamilySize3 True --subset_tag True --nproc int --minFS int --maxFS int --nr_above_bars True/False --output_tabular outptufile_name_tabular --output_pdf outputfile_name_pdf

import argparse
import itertools
import operator
import sys
from collections import Counter
from functools import partial
from multiprocessing.pool import Pool

import matplotlib.pyplot as plt
import numpy
from matplotlib.backends.backend_pdf import PdfPages

plt.switch_backend('agg')


def plotFSDwithHD2(familySizeList1, maximumXFS, minimumXFS, originalCounts,
                   title_file1, subtitle, pdf, relative=False, diff=True):
    if diff is False:
        colors = ["#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4"]
        labels = ["HD=1", "HD=2", "HD=3", "HD=4", "HD=5-8", "HD>8"]
    else:
        colors = ["#93A6AB", "#403C14", "#731E41", "#BAB591", "#085B6F", "#E8AA35", "#726C66"]
        if relative is True:
            labels = ["d=0", "d=0.1", "d=0.2", "d=0.3", "d=0.4", "d=0.5-0.8", "d>0.8"]
        else:
            labels = ["d=0", "d=1", "d=2", "d=3", "d=4", "d=5-8", "d>8"]

    fig = plt.figure(figsize=(6, 7))
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.1)
    p1 = numpy.bincount(numpy.concatenate((familySizeList1)))
    maximumY = numpy.amax(p1)

    if len(range(minimumXFS, maximumXFS)) == 0:
        range1 = range(minimumXFS - 1, minimumXFS + 2)
    else:
        range1 = range(0, maximumXFS + 2)
    counts = plt.hist(familySizeList1, label=labels,
                      color=colors, stacked=True,
                      rwidth=0.8, alpha=1, align="left",
                      edgecolor="None", bins=range1)
    plt.legend(loc='upper right', fontsize=14, frameon=True, bbox_to_anchor=(1.45, 1))

    # plt.title(title_file1, fontsize=12)
    plt.suptitle(subtitle, y=1, x=0.5, fontsize=14)
    plt.xlabel("Family size", fontsize=14)
    plt.ylabel("Absolute Frequency", fontsize=14)

    ticks = numpy.arange(0, maximumXFS + 1, 1)
    ticks1 = map(str, ticks)
    if maximumXFS >= 20:
        ticks1[len(ticks1) - 1] = ">=20"
    plt.xticks(numpy.array(ticks), ticks1)
    [l.set_visible(False) for (i, l) in enumerate(ax.get_xticklabels()) if i % 5 != 0]

    plt.xlim((0, maximumXFS + 1))
    if len(numpy.concatenate(familySizeList1)) != 0:
        plt.ylim((0, max(numpy.bincount(numpy.concatenate(familySizeList1))) * 1.1))

    plt.ylim((0, maximumY * 1.2))
    legend = "\nmax. family size: \nabsolute frequency: \nrelative frequency: "
    plt.text(0.15, -0.08, legend, size=12, transform=plt.gcf().transFigure)

    count = numpy.bincount(originalCounts)  # original counts
    legend1 = "{}\n{}\n{:.5f}".format(max(originalCounts), count[len(count) - 1], float(count[len(count) - 1]) / sum(count))
    plt.text(0.5, -0.08, legend1, size=12, transform=plt.gcf().transFigure)
    legend3 = "singletons\n{:,}\n{:.5f}".format(int(counts[0][len(counts[0]) - 1][1]), float(counts[0][len(counts[0]) - 1][1]) / sum(counts[0][len(counts[0]) - 1]))
    plt.text(0.7, -0.08, legend3, transform=plt.gcf().transFigure, size=12)
    plt.grid(b=True, which='major', color='#424242', linestyle=':')

    pdf.savefig(fig, bbox_inches="tight")
    plt.close("all")


def plotHDwithFSD(list1, maximumX, minimumX, subtitle, lenTags, title_file1, pdf, xlabel, relative=False, nr_above_bars=True):
    if relative is True:
        step = 0.1
    else:
        step = 1

    fig = plt.figure(figsize=(6, 8))
    plt.subplots_adjust(bottom=0.1)
    con_list1 = numpy.concatenate(list1)
    p1 = numpy.array([v for k, v in sorted(Counter(con_list1).iteritems())])
    maximumY = numpy.amax(p1)

    if relative is True:  # relative difference
        bin1 = numpy.arange(-1, maximumX + 0.2, 0.1)
    else:
        bin1 = maximumX + 1

    counts = plt.hist(list1, bins=bin1, edgecolor='black', linewidth=1,
                      label=["FS=1", "FS=2", "FS=3", "FS=4", "FS=5-10",
                             "FS>10"], rwidth=0.8,
                      color=["#808080", "#FFFFCC", "#FFBF00", "#DF0101", "#0431B4", "#86B404"],
                      stacked=True, alpha=1,
                      align="left",
                      range=(0, maximumX + 1))
    plt.legend(loc='upper right', fontsize=14, frameon=True, bbox_to_anchor=(1.45, 1))
    bins = counts[1]  # width of bins
    counts = numpy.array(map(int, counts[0][5]))
    plt.suptitle(subtitle, y=1, x=0.5, fontsize=14)
    # plt.title(title_file1, fontsize=12)
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel("Absolute Frequency", fontsize=14)

    plt.grid(b=True, which='major', color='#424242', linestyle=':')
    plt.axis((minimumX - step, maximumX + step, 0, numpy.amax(counts) + sum(counts) * 0.1))
    plt.xticks(numpy.arange(0, maximumX + step, step))

    plt.ylim((0, maximumY * 1.2))

    if nr_above_bars is True:
        bin_centers = -0.4 * numpy.diff(bins) + bins[:-1]
        for x_label, label in zip(counts, bin_centers):  # labels for values
            if x_label == 0:
                continue
            else:
                plt.annotate("{:,}\n{:.3f}".format(x_label, float(x_label) / sum(counts), 1),
                             xy=(label, x_label + len(con_list1) * 0.01),
                             xycoords="data", color="#000066", fontsize=10)

    legend = "sample size= {:,} against {:,}".format(sum(counts), lenTags)
    plt.text(0.14, -0.01, legend, size=12, transform=plt.gcf().transFigure)

    pdf.savefig(fig, bbox_inches="tight")
    plt.close("all")
    plt.clf()


def plotHDwithinSeq_Sum2(sum1, sum1min, sum2, sum2min, min_value, lenTags, title_file1, pdf):
    fig = plt.figure(figsize=(6, 8))
    plt.subplots_adjust(bottom=0.1)

    ham_partial = [sum1, sum1min, sum2, sum2min, numpy.array(min_value)]  # new hd within tags

    maximumX = numpy.amax(numpy.concatenate(ham_partial))
    minimumX = numpy.amin(numpy.concatenate(ham_partial))
    maximumY = numpy.amax(numpy.array(numpy.concatenate(map(lambda (x): numpy.bincount(x), ham_partial))))

    if len(range(minimumX, maximumX)) == 0:
        range1 = minimumX
    else:
        range1 = range(minimumX, maximumX + 2)

    plt.hist(ham_partial, align="left", rwidth=0.8, stacked=False, label=[ "HD a", "HD b'", "HD b", "HD a'", "HD a+b"], bins=range1, color=["#58ACFA", "#0404B4", "#FE642E", "#B40431", "#585858"], edgecolor='black', linewidth=1)

    plt.legend(loc='upper right', fontsize=14, frameon=True, bbox_to_anchor=(1.55, 1))
    plt.suptitle('Hamming distances within tags', fontsize=14)
    # plt.title(title_file1, fontsize=12)
    plt.xlabel("HD", fontsize=14)
    plt.ylabel("Absolute Frequency", fontsize=14)
    plt.grid(b=True, which='major', color='#424242', linestyle=':')

    plt.axis((minimumX - 1, maximumX + 1, 0, maximumY * 1.2))
    plt.xticks(numpy.arange(0, maximumX + 1, 1.0))
    # plt.ylim(0, maximumY * 1.2)

    legend = "sample size= {:,} against {:,}".format(sum(ham_partial[4]), lenTags)
    plt.text(0.14, -0.01, legend, size=12, transform=plt.gcf().transFigure)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close("all")
    plt.clf()


def createTableFSD2(list1, diff=True):
    selfAB = numpy.concatenate(list1)
    uniqueFS = numpy.unique(selfAB)
    nr = numpy.arange(0, len(uniqueFS), 1)
    if diff is False:
        count = numpy.zeros((len(uniqueFS), 6))
    else:
        count = numpy.zeros((len(uniqueFS), 7))

    state = 1
    for i in list1:
        counts = list(Counter(i).items())
        hd = [item[0] for item in counts]
        c = [item[1] for item in counts]
        table = numpy.column_stack((hd, c))
        if len(table) == 0:
            state = state + 1
            continue
        else:
            if state == 1:
                for i, l in zip(uniqueFS, nr):
                    for j in table:
                        if j[0] == uniqueFS[l]:
                            count[l, 0] = j[1]
            if state == 2:
                for i, l in zip(uniqueFS, nr):
                    for j in table:
                        if j[0] == uniqueFS[l]:
                            count[l, 1] = j[1]

            if state == 3:
                for i, l in zip(uniqueFS, nr):
                    for j in table:
                        if j[0] == uniqueFS[l]:
                            count[l, 2] = j[1]

            if state == 4:
                for i, l in zip(uniqueFS, nr):
                    for j in table:
                        if j[0] == uniqueFS[l]:
                            count[l, 3] = j[1]

            if state == 5:
                for i, l in zip(uniqueFS, nr):
                    for j in table:
                        if j[0] == uniqueFS[l]:
                            count[l, 4] = j[1]

            if state == 6:
                for i, l in zip(uniqueFS, nr):
                    for j in table:
                        if j[0] == uniqueFS[l]:
                            count[l, 5] = j[1]

            if state == 7:
                for i, l in zip(uniqueFS, nr):
                    for j in table:
                        if j[0] == uniqueFS[l]:
                            count[l, 6] = j[1]
            state = state + 1

        sumRow = count.sum(axis=1)
        sumCol = count.sum(axis=0)

    uniqueFS = uniqueFS.astype(str)
    if uniqueFS[len(uniqueFS) - 1] == "20":
        uniqueFS[len(uniqueFS) - 1] = ">20"

    first = ["FS={}".format(i) for i in uniqueFS]
    final = numpy.column_stack((first, count, sumRow))

    return (final, sumCol)


def createFileFSD2(summary, sumCol, overallSum, output_file, name, sep, rel=False, diff=True):
    output_file.write(name)
    output_file.write("\n")
    if diff is False:
        output_file.write("{}HD=1{}HD=2{}HD=3{}HD=4{}HD=5-8{}HD>8{}sum{}\n".format(sep, sep, sep, sep, sep, sep, sep, sep))
    else:
        if rel is False:
            output_file.write("{}diff=0{}diff=1{}diff=2{}diff=3{}diff=4{}diff=5-8{}diff>8{}sum{}\n".format(sep, sep, sep, sep, sep, sep, sep, sep, sep))
        else:
            output_file.write("{}diff=0{}diff=0.1{}diff=0.2{}diff=0.3{}diff=0.4{}diff=0.5-0.8{}diff>0.8{}sum{}\n".format(sep, sep, sep, sep, sep, sep, sep, sep, sep))

    for item in summary:
        for nr in item:
            if "FS" not in nr and "diff" not in nr:
                nr = nr.astype(float)
                nr = nr.astype(int)
            output_file.write("{}{}".format(nr, sep))
        output_file.write("\n")
    output_file.write("sum{}".format(sep))
    sumCol = map(int, sumCol)
    for el in sumCol:
        output_file.write("{}{}".format(el, sep))
    output_file.write("{}{}".format(overallSum.astype(int), sep))
    output_file.write("\n\n")


def createTableHD(list1, row_label):
    selfAB = numpy.concatenate(list1)
    uniqueHD = numpy.unique(selfAB)
    nr = numpy.arange(0, len(uniqueHD), 1)
    count = numpy.zeros((len(uniqueHD), 6))
    state = 1
    for i in list1:
        counts = list(Counter(i).items())
        hd = [item[0] for item in counts]
        c = [item[1] for item in counts]
        table = numpy.column_stack((hd, c))
        if len(table) == 0:
            state = state + 1
            continue
        else:
            if state == 1:
                for i, l in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 0] = j[1]
            if state == 2:
                for i, l in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 1] = j[1]

            if state == 3:
                for i, l in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 2] = j[1]

            if state == 4:
                for i, l in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 3] = j[1]

            if state == 5:
                for i, l in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 4] = j[1]

            if state == 6:
                for i, l in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 5] = j[1]
            state = state + 1

        sumRow = count.sum(axis=1)
        sumCol = count.sum(axis=0)
        first = ["{}{}".format(row_label, i) for i in uniqueHD]
        final = numpy.column_stack((first, count, sumRow))

    return (final, sumCol)


def createTableHDwithTags(list1):
    selfAB = numpy.concatenate(list1)
    uniqueHD = numpy.unique(selfAB)
    nr = numpy.arange(0, len(uniqueHD), 1)
    count = numpy.zeros((len(uniqueHD), 5))

    state = 1
    for i in list1:
        counts = list(Counter(i).items())
        hd = [item[0] for item in counts]
        c = [item[1] for item in counts]
        table = numpy.column_stack((hd, c))
        if len(table) == 0:
            state = state + 1
            continue
        else:
            if state == 1:
                for i, l in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 0] = j[1]
            if state == 2:
                for i, l in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 1] = j[1]
            if state == 3:
                for i, l in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 2] = j[1]
            if state == 4:
                for i, l in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 3] = j[1]
            if state == 5:
                for i, l in zip(uniqueHD, nr):
                    for j in table:
                        if j[0] == uniqueHD[l]:
                            count[l, 4] = j[1]

            state = state + 1

        sumRow = count.sum(axis=1)
        sumCol = count.sum(axis=0)
        first = ["HD={}".format(i) for i in uniqueHD]
        final = numpy.column_stack((first, count, sumRow))

    return (final, sumCol)


def createFileHD(summary, sumCol, overallSum, output_file, name, sep):
    output_file.write(name)
    output_file.write("\n")
    output_file.write("{}FS=1{}FS=2{}FS=3{}FS=4{}FS=5-10{}FS>10{}sum{}\n".format(sep, sep, sep, sep, sep, sep, sep, sep))
    for item in summary:
        for nr in item:
            if "HD" not in nr and "diff" not in nr:
                nr = nr.astype(float)
                nr = nr.astype(int)
            output_file.write("{}{}".format(nr, sep))
        output_file.write("\n")
    output_file.write("sum{}".format(sep))
    sumCol = map(int, sumCol)
    for el in sumCol:
        output_file.write("{}{}".format(el, sep))
    output_file.write("{}{}".format(overallSum.astype(int), sep))
    output_file.write("\n\n")


def createFileHDwithinTag(summary, sumCol, overallSum, output_file, name, sep):
    output_file.write(name)
    output_file.write("\n")
    output_file.write("{}HD a{}HD b'{}HD b{}HD a'{}HD a+b{}sum{}\n".format(sep, sep, sep, sep, sep, sep, sep))
    for item in summary:
        for nr in item:
            if "HD" not in nr:
                nr = nr.astype(float)
                nr = nr.astype(int)
            output_file.write("{}{}".format(nr, sep))
        output_file.write("\n")
    output_file.write("sum{}".format(sep))
    sumCol = map(int, sumCol)
    for el in sumCol:
        output_file.write("{}{}".format(el, sep))
    output_file.write("{}{}".format(overallSum.astype(int), sep))
    output_file.write("\n\n")


def hamming(array1, array2):
    res = 99 * numpy.ones(len(array1))
    i = 0
    array2 = numpy.unique(array2)  # remove duplicate sequences to decrease running time
    for a in array1:
        dist = numpy.array([sum(itertools.imap(operator.ne, a, b)) for b in array2])  # fastest
        res[i] = numpy.amin(dist[dist > 0])  # pick min distance greater than zero
        # print(i)
        i += 1
    return res


def hamming_difference(array1, array2, mate_b):
    array2 = numpy.unique(array2)  # remove duplicate sequences to decrease running time
    array1_half = numpy.array([i[0:(len(i)) / 2] for i in array1])  # mate1 part1
    array1_half2 = numpy.array([i[len(i) / 2:len(i)] for i in array1])  # mate1 part 2

    array2_half = numpy.array([i[0:(len(i)) / 2] for i in array2])  # mate2 part1
    array2_half2 = numpy.array([i[len(i) / 2:len(i)] for i in array2])  # mate2 part2

    # diff11 = 999 * numpy.ones(len(array2))
    # relativeDiffList = 999 * numpy.ones(len(array2))
    # ham1 = 999 * numpy.ones(len(array2))
    # ham2 = 999 * numpy.ones(len(array2))
    # min_valueList = 999 * numpy.ones(len(array2))
    # min_tagsList = 999 * numpy.ones(len(array2))
    # diff11_zeros = 999 * numpy.ones(len(array2))
    # min_tagsList_zeros = 999 * numpy.ones(len(array2))

    diff11 = []
    relativeDiffList = []
    ham1 = []
    ham2 = []
    ham1min = []
    ham2min = []
    min_valueList = []
    min_tagsList = []
    diff11_zeros = []
    min_tagsList_zeros = []
    i = 0  # counter, only used to see how many HDs of tags were already calculated
    if mate_b is False:  # HD calculation for all a's
        half1_mate1 = array1_half
        half2_mate1 = array1_half2
        half1_mate2 = array2_half
        half2_mate2 = array2_half2
    elif mate_b is True:  # HD calculation for all b's
        half1_mate1 = array1_half2
        half2_mate1 = array1_half
        half1_mate2 = array2_half2
        half2_mate2 = array2_half

    for a, b, tag in zip(half1_mate1, half2_mate1, array1):
        # exclude identical tag from array2, to prevent comparison to itself
        sameTag = numpy.where(array2 == tag)
        indexArray2 = numpy.arange(0, len(array2), 1)
        index_withoutSame = numpy.delete(indexArray2, sameTag)  # delete identical tag from the data

        # all tags without identical tag
        array2_half_withoutSame = half1_mate2[index_withoutSame]
        array2_half2_withoutSame = half2_mate2[index_withoutSame]
        # array2_withoutSame = array2[index_withoutSame]  # whole tag (=not splitted into 2 halfs)

        dist = numpy.array([sum(itertools.imap(operator.ne, a, c)) for c in
                            array2_half_withoutSame])  # calculate HD of "a" in the tag to all "a's" or "b" in the tag to all "b's"
        min_index = numpy.where(dist == dist.min())  # get index of min HD
        min_value = dist[min_index]  # get minimum HDs
        min_tag_half2 = array2_half2_withoutSame[min_index]  # get all "b's" of the tag or all "a's" of the tag with minimum HD
        # min_tag = array2_withoutSame[min_index]  # get whole tag with min HD

        dist2 = numpy.array([sum(itertools.imap(operator.ne, b, e)) for e in
                             min_tag_half2])  # calculate HD of "b" to all "b's" or "a" to all "a's"
        for d, d2 in zip(min_value, dist2):
            if mate_b is True:  # half2, corrects the variable of the HD from both halfs if it is a or b
                ham2.append(d)
                ham2min.append(d2)
            else:  # half1, corrects the variable of the HD from both halfs if it is a or b
                ham1.append(d)
                ham1min.append(d2)

            min_valueList.append(d + d2)
            min_tagsList.append(tag)
            difference1 = abs(d - d2)
            diff11.append(difference1)
            rel_difference = round(float(difference1) / (d + d2), 1)
            relativeDiffList.append(rel_difference)

            # tags which have identical parts:
            if d == 0 or d2 == 0:
                min_tagsList_zeros.append(tag)
                difference1_zeros = abs(d - d2)
                diff11_zeros.append(difference1_zeros)
        i += 1

    # print(i)
    # diff11 = [st for st in diff11 if st != 999]
    # ham1 = [st for st in ham1 if st != 999]
    # ham2 = [st for st in ham2 if st != 999]
    # min_valueList = [st for st in min_valueList if st != 999]
    # min_tagsList = [st for st in min_tagsList if st != 999]
    # relativeDiffList = [st for st in relativeDiffList if st != 999]
    # diff11_zeros = [st for st in diff11_zeros if st != 999]
    # min_tagsList_zeros = [st for st in min_tagsList_zeros if st != 999]

    return ([diff11, ham1, ham2, min_valueList, min_tagsList, relativeDiffList, diff11_zeros, min_tagsList_zeros, ham1min, ham2min])


def readFileReferenceFree(file):
    with open(file, 'r') as dest_f:
        data_array = numpy.genfromtxt(dest_f, skip_header=0, delimiter='\t', comments='#', dtype='string')
        integers = numpy.array(data_array[:, 0]).astype(int)
        return(integers, data_array)


def hammingDistanceWithFS(fs, ham):
    fs = numpy.asarray(fs)
    maximum = max(ham)
    minimum = min(ham)
    ham = numpy.asarray(ham)

    singletons = numpy.where(fs == 1)[0]
    data = ham[singletons]

    hd2 = numpy.where(fs == 2)[0]
    data2 = ham[hd2]

    hd3 = numpy.where(fs == 3)[0]
    data3 = ham[hd3]

    hd4 = numpy.where(fs == 4)[0]
    data4 = ham[hd4]

    hd5 = numpy.where((fs >= 5) & (fs <= 10))[0]
    data5 = ham[hd5]

    hd6 = numpy.where(fs > 10)[0]
    data6 = ham[hd6]

    list1 = [data, data2, data3, data4, data5, data6]
    return(list1, maximum, minimum)


def familySizeDistributionWithHD(fs, ham, diff=False, rel=True):
    hammingDistances = numpy.unique(ham)
    fs = numpy.asarray(fs)

    ham = numpy.asarray(ham)
    bigFamilies2 = numpy.where(fs > 19)[0]
    if len(bigFamilies2) != 0:
        fs[bigFamilies2] = 20
    maximum = max(fs)
    minimum = min(fs)
    if diff is True:
        hd0 = numpy.where(ham == 0)[0]
        data0 = fs[hd0]

    if rel is True:
        hd1 = numpy.where(ham == 0.1)[0]
    else:
        hd1 = numpy.where(ham == 1)[0]
    data = fs[hd1]

    if rel is True:
        hd2 = numpy.where(ham == 0.2)[0]
    else:
        hd2 = numpy.where(ham == 2)[0]
    data2 = fs[hd2]

    if rel is True:
        hd3 = numpy.where(ham == 0.3)[0]
    else:
        hd3 = numpy.where(ham == 3)[0]
    data3 = fs[hd3]

    if rel is True:
        hd4 = numpy.where(ham == 0.4)[0]
    else:
        hd4 = numpy.where(ham == 4)[0]
    data4 = fs[hd4]

    if rel is True:
        hd5 = numpy.where((ham >= 0.5) & (ham <= 0.8))[0]
    else:
        hd5 = numpy.where((ham >= 5) & (ham <= 8))[0]
    data5 = fs[hd5]

    if rel is True:
        hd6 = numpy.where(ham > 0.8)[0]
    else:
        hd6 = numpy.where(ham > 8)[0]
    data6 = fs[hd6]

    if diff is True:
        list1 = [data0, data, data2, data3, data4, data5, data6]
    else:
        list1 = [data, data2, data3, data4, data5, data6]

    return(list1, hammingDistances, maximum, minimum)


def make_argparser():
    parser = argparse.ArgumentParser(description='Hamming distance analysis of duplex sequencing data')
    parser.add_argument('--inputFile',
                        help='Tabular File with three columns: ab or ba, tag and family size.')
    parser.add_argument('--inputName1')
    parser.add_argument('--inputFile2', default=None,
                        help='Tabular File with three columns: ab or ba, tag and family size.')
    parser.add_argument('--inputName2')
    parser.add_argument('--sample_size', default=1000, type=int,
                        help='Sample size of Hamming distance analysis.')
    parser.add_argument('--subset_tag', default=0, type=int,
                        help='The tag is shortened to the given number.')
    parser.add_argument('--nproc', default=4, type=int,
                        help='The tool runs with the given number of processors.')
    parser.add_argument('--only_DCS', action="store_false",
                        help='Only tags of the DCSs are included in the HD analysis')

    parser.add_argument('--minFS', default=1, type=int,
                        help='Only tags, which have a family size greater or equal than specified, are included in the HD analysis')
    parser.add_argument('--maxFS', default=0, type=int,
                        help='Only tags, which have a family size smaller or equal than specified, are included in the HD analysis')
    parser.add_argument('--nr_above_bars', action="store_true",
                        help='If no, values above bars in the histrograms are removed')

    parser.add_argument('--output_tabular', default="data.tabular", type=str,
                        help='Name of the tabular file.')
    parser.add_argument('--output_pdf', default="data.pdf", type=str,
                        help='Name of the pdf file.')
    parser.add_argument('--output_pdf2', default="data2.pdf", type=str,
                        help='Name of the pdf file.')
    parser.add_argument('--output_tabular2', default="data2.tabular", type=str,
                        help='Name of the tabular file.')

    return parser


def Hamming_Distance_Analysis(argv):
    parser = make_argparser()
    args = parser.parse_args(argv[1:])

    file1 = args.inputFile
    name1 = args.inputName1

    file2 = args.inputFile2
    name2 = args.inputName2

    index_size = args.sample_size
    title_savedFile_pdf = args.output_pdf
    title_savedFile_pdf2 = args.output_pdf2

    title_savedFile_csv = args.output_tabular
    title_savedFile_csv2 = args.output_tabular2

    sep = "\t"
    onlyDuplicates = args.only_DCS
    minFS = args.minFS
    maxFS = args.maxFS
    nr_above_bars = args.nr_above_bars

    subset = args.subset_tag
    nproc = args.nproc

    # input checks
    if index_size < 0:
        print("index_size is a negative integer.")
        exit(2)

    if nproc <= 0:
        print("nproc is smaller or equal zero")
        exit(3)

    if subset < 0:
        print("subset_tag is smaller or equal zero.")
        exit(5)

    # PLOT
    plt.rcParams['axes.facecolor'] = "E0E0E0"  # grey background color
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14
    plt.rcParams['patch.edgecolor'] = "#000000"
    plt.rc('figure', figsize=(11.69, 8.27))  # A4 format

    if file2 != str(None):
        files = [file1, file2]
        name1 = name1.split(".tabular")[0]
        name2 = name2.split(".tabular")[0]
        names = [name1, name2]
        pdf_files = [title_savedFile_pdf, title_savedFile_pdf2]
        csv_files = [title_savedFile_csv, title_savedFile_csv2]
    else:
        files = [file1]
        name1 = name1.split(".tabular")[0]
        names = [name1]
        pdf_files = [title_savedFile_pdf]
        csv_files = [title_savedFile_csv]

    for f, name_file, pdf_f, csv_f in zip(files, names, pdf_files, csv_files):
        with open(csv_f, "w") as output_file, PdfPages(pdf_f) as pdf:
            print("dataset: ", name_file)
            integers, data_array = readFileReferenceFree(f)
            data_array = numpy.array(data_array)
            int_f = numpy.array(data_array[:, 0]).astype(int)
            data_array = data_array[numpy.where(int_f >= minFS)]
            integers = integers[integers >= minFS]

            # select family size for tags
            if maxFS > 0:
                int_f2 = numpy.array(data_array[:, 0]).astype(int)
                data_array = data_array[numpy.where(int_f2 <= maxFS)]
                integers = integers[integers <= maxFS]

            print("min FS", min(integers))
            print("max FS", max(integers))

            tags = data_array[:, 2]
            seq = data_array[:, 1]

            if onlyDuplicates is True:
                # find all unique tags and get the indices for ALL tags, but only once
                u, index_unique, c = numpy.unique(numpy.array(seq), return_counts=True, return_index=True)
                d = u[c > 1]

                # get family sizes, tag for duplicates
                duplTags_double = integers[numpy.in1d(seq, d)]
                duplTags = duplTags_double[0::2]  # ab of DCS
                duplTagsBA = duplTags_double[1::2]  # ba of DCS

                duplTags_tag = tags[numpy.in1d(seq, d)][0::2]  # ab
                duplTags_seq = seq[numpy.in1d(seq, d)][0::2]  # ab - tags

                data_array = numpy.column_stack((duplTags, duplTags_seq))
                data_array = numpy.column_stack((data_array, duplTags_tag))
                integers = numpy.array(data_array[:, 0]).astype(int)
                print("DCS in whole dataset", len(data_array))

            # HD analysis for a subset of the tag
            if subset > 0:
                tag1 = numpy.array([i[0:(len(i)) / 2] for i in data_array[:, 1]])
                tag2 = numpy.array([i[len(i) / 2:len(i)] for i in data_array[:, 1]])

                flanking_region_float = float((len(tag1[0]) - subset)) / 2
                flanking_region = int(flanking_region_float)
                if flanking_region_float % 2 == 0:
                    tag1_shorten = numpy.array([i[flanking_region:len(i) - flanking_region] for i in tag1])
                    tag2_shorten = numpy.array([i[flanking_region:len(i) - flanking_region] for i in tag2])
                else:
                    flanking_region_rounded = int(round(flanking_region, 1))
                    flanking_region_rounded_end = len(tag1[0]) - subset - flanking_region_rounded
                    tag1_shorten = numpy.array(
                        [i[flanking_region:len(i) - flanking_region_rounded_end] for i in tag1])
                    tag2_shorten = numpy.array(
                        [i[flanking_region:len(i) - flanking_region_rounded_end] for i in tag2])

                data_array_tag = numpy.array([i + j for i, j in zip(tag1_shorten, tag2_shorten)])
                data_array = numpy.column_stack((data_array[:, 0], data_array_tag, data_array[:, 2]))

            print("length of tag= ", len(data_array[0, 1]))
            # select sample: if no size given --> all vs. all comparison
            if index_size == 0:
                result = numpy.arange(0, len(data_array), 1)
            else:
                result = numpy.random.choice(len(integers), size=index_size, replace=False)  # array of random sequences of size=index.size

            # with open("index_result1_{}.pkl".format(app_f), "wb") as o:
            #     pickle.dump(result, o, pickle.HIGHEST_PROTOCOL)

            # comparison random tags to whole dataset
            result1 = data_array[result, 1]  # random tags
            result2 = data_array[:, 1]  # all tags
            print("size of the whole dataset= ", len(result2))
            print("sample size= ", len(result1))

            # HD analysis of whole tag
            proc_pool = Pool(nproc)
            chunks_sample = numpy.array_split(result1, nproc)
            ham = proc_pool.map(partial(hamming, array2=result2), chunks_sample)
            proc_pool.close()
            proc_pool.join()
            ham = numpy.concatenate(ham).astype(int)
            # with open("HD_whole dataset_{}.txt".format(app_f), "w") as output_file1:
            # for h, tag in zip(ham, result1):
            #     output_file1.write("{}\t{}\n".format(tag, h))

            # HD analysis for chimeric reads
            proc_pool_b = Pool(nproc)
            diff_list_a = proc_pool_b.map(partial(hamming_difference, array2=result2, mate_b=False), chunks_sample)
            diff_list_b = proc_pool_b.map(partial(hamming_difference, array2=result2, mate_b=True), chunks_sample)
            proc_pool_b.close()
            proc_pool_b.join()
            diff = numpy.concatenate((numpy.concatenate([item[0] for item in diff_list_a]),
                                      numpy.concatenate([item_b[0] for item_b in diff_list_b]))).astype(int)
            HDhalf1 = numpy.concatenate((numpy.concatenate([item[1] for item in diff_list_a]),
                                         numpy.concatenate([item_b[1] for item_b in diff_list_b]))).astype(int)
            HDhalf2 = numpy.concatenate((numpy.concatenate([item[2] for item in diff_list_a]),
                                         numpy.concatenate([item_b[2] for item_b in diff_list_b]))).astype(int)
            minHDs = numpy.concatenate((numpy.concatenate([item[3] for item in diff_list_a]),
                                        numpy.concatenate([item_b[3] for item_b in diff_list_b]))).astype(int)
            minHD_tags = numpy.concatenate((numpy.concatenate([item[4] for item in diff_list_a]),
                                            numpy.concatenate([item_b[4] for item_b in diff_list_b])))
            rel_Diff = numpy.concatenate((numpy.concatenate([item[5] for item in diff_list_a]),
                                          numpy.concatenate([item_b[5] for item_b in diff_list_b])))
            diff_zeros = numpy.concatenate((numpy.concatenate([item[6] for item in diff_list_a]),
                                            numpy.concatenate([item_b[6] for item_b in diff_list_b]))).astype(int)
            minHD_tags_zeros = numpy.concatenate((numpy.concatenate([item[7] for item in diff_list_a]),
                                                  numpy.concatenate([item_b[7] for item_b in diff_list_b])))
            HDhalf1min = numpy.concatenate((numpy.concatenate([item[8] for item in diff_list_a]), numpy.concatenate([item_b[8] for item_b in diff_list_b]))).astype(int)
            HDhalf2min = numpy.concatenate((numpy.concatenate([item[9] for item in diff_list_a]),
                                            numpy.concatenate([item_b[9] for item_b in diff_list_b]))).astype(int)

            lenTags = len(data_array)

            quant = numpy.array(data_array[result, 0]).astype(int)  # family size for sample of tags
            seq = numpy.array(data_array[result, 1])  # tags of sample
            ham = numpy.asarray(ham)  # HD for sample of tags

            if onlyDuplicates is True:  # ab and ba strands of DCSs
                quant = numpy.concatenate((quant, duplTagsBA[result]))
                seq = numpy.tile(seq, 2)
                ham = numpy.tile(ham, 2)

            # prepare data for different kinds of plots
            # distribution of FSs separated after HD
            familySizeList1, hammingDistances, maximumXFS, minimumXFS = familySizeDistributionWithHD(quant, ham, rel=False)
            list1, maximumX, minimumX = hammingDistanceWithFS(quant, ham)  # histogram of HDs separated after FS

            # get FS for all tags with min HD of analysis of chimeric reads
            # there are more tags than sample size in the plot, because one tag can have multiple minimas
            seqDic = dict(zip(seq, quant))
            lst_minHD_tags = []
            for i in minHD_tags:
                lst_minHD_tags.append(seqDic.get(i))

            # histogram with absolute and relative difference between HDs of both parts of the tag
            listDifference1, maximumXDifference, minimumXDifference = hammingDistanceWithFS(lst_minHD_tags, diff)
            listRelDifference1, maximumXRelDifference, minimumXRelDifference = hammingDistanceWithFS(lst_minHD_tags,
                                                                                                     rel_Diff)

            # family size distribution separated after the difference between HDs of both parts of the tag
            familySizeList1_diff, hammingDistances_diff, maximumXFS_diff, minimumXFS_diff = familySizeDistributionWithHD(
                lst_minHD_tags, diff, diff=True, rel=False)
            familySizeList1_reldiff, hammingDistances_reldiff, maximumXFS_reldiff, minimumXFS_reldiff = familySizeDistributionWithHD(
                lst_minHD_tags, rel_Diff, diff=True, rel=True)

            # chimeric read analysis: tags which have HD=0 in one of the halfs
            if len(minHD_tags_zeros) != 0:
                lst_minHD_tags_zeros = []
                for i in minHD_tags_zeros:
                    lst_minHD_tags_zeros.append(seqDic.get(i))  # get family size for tags of chimeric reads

            # histogram with HD of non-identical half
            listDifference1_zeros, maximumXDifference_zeros, minimumXDifference_zeros = hammingDistanceWithFS(lst_minHD_tags_zeros, diff_zeros)
            # family size distribution of non-identical half
            familySizeList1_diff_zeros, hammingDistances_diff_zeros, maximumXFS_diff_zeros, minimumXFS_diff_zeros = familySizeDistributionWithHD(lst_minHD_tags_zeros, diff_zeros, diff=False, rel=False)

            # plot Hamming Distance with Family size distribution
            plotHDwithFSD(list1=list1, maximumX=maximumX, minimumX=minimumX, pdf=pdf, subtitle="Hamming distance separated by family size", title_file1=name_file, lenTags=lenTags, xlabel="HD", nr_above_bars=nr_above_bars)

            # Plot FSD with separation after
            plotFSDwithHD2(familySizeList1, maximumXFS, minimumXFS,
                           originalCounts=quant, subtitle="Family size distribution separated by Hamming distance",
                           pdf=pdf, relative=False, title_file1=name_file, diff=False)

            # Plot HD within tags
            plotHDwithinSeq_Sum2(HDhalf1, HDhalf1min, HDhalf2, HDhalf2min, minHDs, pdf=pdf, lenTags=lenTags, title_file1=name_file)

            # Plot difference between HD's separated after FSD
            plotHDwithFSD(listDifference1, maximumXDifference, minimumXDifference, pdf=pdf,
                          subtitle="Delta Hamming distance within tags",
                          title_file1=name_file, lenTags=lenTags,
                          xlabel="absolute delta HD", relative=False, nr_above_bars=nr_above_bars)

            plotHDwithFSD(listRelDifference1, maximumXRelDifference, minimumXRelDifference, pdf=pdf,
                          subtitle="Chimera Analysis: relative delta Hamming distances",
                          title_file1=name_file, lenTags=lenTags,
                          xlabel="relative delta HD", relative=True, nr_above_bars=nr_above_bars)

            # plots for chimeric reads
            if len(minHD_tags_zeros) != 0:
                # HD
                plotHDwithFSD(listDifference1_zeros, maximumXDifference_zeros, minimumXDifference_zeros, pdf=pdf,
                              subtitle="Hamming distance of the non-identical half of chimeras",
                              title_file1=name_file, lenTags=lenTags, xlabel="HD", relative=False, nr_above_bars=nr_above_bars)

            # print all data to a CSV file
            # HD
            summary, sumCol = createTableHD(list1, "HD=")
            overallSum = sum(sumCol)  # sum of columns in table

            # FSD
            summary5, sumCol5 = createTableFSD2(familySizeList1, diff=False)
            overallSum5 = sum(sumCol5)

            # HD of both parts of the tag
            summary9, sumCol9 = createTableHDwithTags([HDhalf1, HDhalf1min, HDhalf2, HDhalf2min, numpy.array(minHDs)])
            overallSum9 = sum(sumCol9)

            # HD
            # absolute difference
            summary11, sumCol11 = createTableHD(listDifference1, "diff=")
            overallSum11 = sum(sumCol11)
            # relative difference and all tags
            summary13, sumCol13 = createTableHD(listRelDifference1, "diff=")
            overallSum13 = sum(sumCol13)

            # chimeric reads
            if len(minHD_tags_zeros) != 0:
                # absolute difference and tags where at least one half has HD=0
                summary15, sumCol15 = createTableHD(listDifference1_zeros, "HD=")
                overallSum15 = sum(sumCol15)

            output_file.write("{}\n".format(name_file))
            output_file.write("number of tags per file{}{:,} (from {:,}) against {:,}\n\n".format(sep, len(
                numpy.concatenate(list1)), lenTags, lenTags))

            # HD
            createFileHD(summary, sumCol, overallSum, output_file,
                         "Hamming distance separated by family size", sep)
            # FSD
            createFileFSD2(summary5, sumCol5, overallSum5, output_file,
                           "Family size distribution separated by Hamming distance", sep,
                           diff=False)

            count = numpy.bincount(quant)
            # output_file.write("{}{}\n".format(sep, name_file))
            output_file.write("\n")
            output_file.write("max. family size:{}{}\n".format(sep, max(quant)))
            output_file.write("absolute frequency:{}{}\n".format(sep, count[len(count) - 1]))
            output_file.write(
                "relative frequency:{}{}\n\n".format(sep, float(count[len(count) - 1]) / sum(count)))

            # HD within tags
            output_file.write(
                "The hamming distances were calculated by comparing each half of all tags against the tag(s) with the minimum Hamming distance per half.\n"
                "It is possible that one tag can have the minimum HD from multiple tags, so the sample size in this calculation differs from the sample size entered by the user.\n")
            output_file.write(
                "actual number of tags with min HD = {:,} (sample size by user = {:,})\n".format(
                    len(numpy.concatenate(listDifference1)), len(numpy.concatenate(list1))))
            output_file.write("length of one part of the tag = {}\n\n".format(len(data_array[0, 1]) / 2))

            createFileHDwithinTag(summary9, sumCol9, overallSum9, output_file,
                                  "Hamming distance of each half in the tag", sep)
            createFileHD(summary11, sumCol11, overallSum11, output_file,
                         "Absolute delta Hamming distances within the tag", sep)
            createFileHD(summary13, sumCol13, overallSum13, output_file,
                         "Chimera analysis: relative delta Hamming distances", sep)

            if len(minHD_tags_zeros) != 0:
                output_file.write(
                    "Chimeras:\nAll tags were filtered: only those tags where at least one half is identical with the half of the min. tag are kept.\nSo the hamming distance of the non-identical half is compared.\n")
                createFileHD(summary15, sumCol15, overallSum15, output_file,
                             "Hamming distances of non-zero half", sep)
            output_file.write("\n")


if __name__ == '__main__':
    sys.exit(Hamming_Distance_Analysis(sys.argv))
