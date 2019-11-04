#!/usr/bin/env python

"""read2mut.py

Author -- Gundula Povysil
Contact -- povysil@bioinf.jku.at

Looks for reads with mutation at known
positions and calculates frequencies and stats.

=======  ==========  =================  ================================
Version  Date        Author             Description
0.2.1    2019-10-27  Gundula Povysil    -
=======  ==========  =================  ================================


USAGE: python read2mut.py --mutFile DCS_Mutations.tabular --bamFile Interesting_Reads.trim.bam
                          --inputJson tag_count_dict.json --sscsJson SSCS_counts.json
                          --outputFile mutant_reads_summary_short_trim.xlsx --thresh 10 --phred 20

"""

import argparse
import sys
import os

from __future__ import division
import numpy as np
import re
import json
import pysam
import itertools
import operator
import xlsxwriter


def make_argparser():
    parser = argparse.ArgumentParser(description='Takes a tabular file with mutations, a BAM file and JSON files as input and prints stats about variants to a user specified output file.')
    parser.add_argument('--mutFile',
                        help='TABULAR file with DCS mutations.')
    parser.add_argument('--bamFile',
                        help='BAM file with aligned raw reads of selected tags (FASTQ created by mut2read.py - trimming with Trimmomatic - alignment with bwa).')
    parser.add_argument('--inputJson',
                        help='JSON file with data collected by mut2read.py.')
    parser.add_argument('--sscsJson',
                        help='JSON file with SSCS counts collected by mut2sscs.py.')
    parser.add_argument('--outputFile',
                        help='Output xlsx file of mutation details.')
    parser.add_argument('--thresh', type=int, default=0,
                        help='Integer threshold for displaying mutations. Only mutations occuring less than thresh times are displayed. Default of 0 displays all.')
    parser.add_argument('--phred', type=int, default=20,
                        help='Integer threshold for Phred score. Only reads higher than this threshold are considered. Default 20.')
    return parser


def read2mut(argv):
    parser = make_argparser()
    args = parser.parse_args(argv[1:])
    file1 = args.mutFile
    file2 = args.bamFile
    json_file = args.inputJson
    sscs_json = args.sscsJson
    outfile = args.outputFile
    thresh = args.thresh
    phred_score = args.phred

    if os.path.isfile(file1) is False:
        print("Error: Could not find '{}'".format(file1))
        exit(0)
    if os.path.isfile(file2) is False:
        print("Error: Could not find '{}'".format(file2))
        exit(1)
    if os.path.isfile(json_file) is False:
        print("Error: Could not find '{}'".format(json_file))
        exit(2)
    if thresh < 0:
        print("Error: Thresh is '{}', but only non-negative integers allowed".format(thresh))
        exit(3)
    if phred_score < 0:
        print("Error: Phred is '{}', but only non-negative integers allowed".format(phred_score))
        exit(4)

    # 1. read mut file
    with open(file1, 'r') as mut:
        mut_array = np.genfromtxt(mut, skip_header=1, delimiter='\t', comments='#', dtype='string')

    # 2. load dicts
    with open(json_file, "r") as f:
        (tag_dict, cvrg_dict) = json.load(f)

    with open(sscs_json, "r") as f:
        (mut_pos_dict, ref_pos_dict) = json.load(f)

    # 3. read bam file
    pysam.index(file2)
    bam = pysam.AlignmentFile(file2, "rb")

    # 4. create mut_dict
    mut_dict = {}
    mut_read_pos_dict = {}
    mut_read_dict = {}
    reads_dict = {}
    if len(mut_array) == 13: 
        mut_array = mut_array.reshape((1, len(mut_array)))

    for m in range(0, len(mut_array[:, 0])):
        print(str(m + 1) + " of " + str(len(mut_array[:, 0])))
        #    for m in range(0, 5):
        chrom = mut_array[m, 1]
        stop_pos = mut_array[m, 2].astype(int)
        chrom_stop_pos = str(chrom) + "#" + str(stop_pos)
        ref = mut_array[m, 9]
        alt = mut_array[m, 10]
        mut_dict[chrom_stop_pos] = {}
        mut_read_pos_dict[chrom_stop_pos] = {}
        reads_dict[chrom_stop_pos] = {}

        for pileupcolumn in bam.pileup(chrom.tobytes(), stop_pos - 2, stop_pos, max_depth=1000000000):
            if pileupcolumn.reference_pos == stop_pos - 1:
                count_alt = 0
                count_ref = 0
                count_indel = 0
                count_n = 0
                count_other = 0
                count_lowq = 0
                n = 0
                print("unfiltered reads=", pileupcolumn.n, "filtered reads=", len(pileupcolumn.pileups),
                      "difference= ", len(pileupcolumn.pileups)-pileupcolumn.n)
                for pileupread in pileupcolumn.pileups:
                    n += 1
                    if not pileupread.is_del and not pileupread.is_refskip:
                        tag = pileupread.alignment.query_name
                        nuc = pileupread.alignment.query_sequence[pileupread.query_position]
                        phred = ord(pileupread.alignment.qual[pileupread.query_position]) - 33
                        if phred < phred_score:
                            nuc = "lowQ"
                        if tag not in mut_dict[chrom_stop_pos]:
                            mut_dict[chrom_stop_pos][tag] = {}
                        if nuc in mut_dict[chrom_stop_pos][tag]:
                            mut_dict[chrom_stop_pos][tag][nuc] += 1
                        else:
                            mut_dict[chrom_stop_pos][tag][nuc] = 1
                        if tag not in mut_read_pos_dict[chrom_stop_pos]:
                            mut_read_pos_dict[chrom_stop_pos][tag] = np.array(pileupread.query_position)
                            reads_dict[chrom_stop_pos][tag] = len(pileupread.alignment.query_sequence)
                        else:
                            mut_read_pos_dict[chrom_stop_pos][tag] = np.append(
                                mut_read_pos_dict[chrom_stop_pos][tag], pileupread.query_position)
                            reads_dict[chrom_stop_pos][tag] = np.append(
                                reads_dict[chrom_stop_pos][tag], len(pileupread.alignment.query_sequence))

                        if nuc == alt:
                            count_alt += 1
                            if tag not in mut_read_dict:
                                mut_read_dict[tag] = {}
                                mut_read_dict[tag][chrom_stop_pos] = alt
                            else:
                                mut_read_dict[tag][chrom_stop_pos] = alt
                        elif nuc == ref:
                            count_ref += 1
                        elif nuc == "N":
                            count_n += 1
                        elif nuc == "lowQ":
                            count_lowq += 1
                        else:
                            count_other += 1
                    else:
                        count_indel += 1

        print("coverage at pos %s = %s, ref = %s, alt = %s, other bases = %s, N = %s, indel = %s, low quality = %s\n"
               % (pileupcolumn.pos, count_ref + count_alt, count_ref, count_alt, count_other, count_n, count_indel,
                  count_lowq))

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            pure_tag = read.query_name[:-5]
            nuc = "na"
            for key in tag_dict[pure_tag].keys():
                if key not in mut_dict:
                    mut_dict[key] = {}
                if read.query_name not in mut_dict[key]:
                    mut_dict[key][read.query_name] = {}
                if nuc in mut_dict[key][read.query_name]:
                    mut_dict[key][read.query_name][nuc] += 1
                else:
                    mut_dict[key][read.query_name][nuc] = 1
    bam.close()

    # 5. create pure_tags_dict
    pure_tags_dict = {}
    for key1, value1 in sorted(mut_dict.items()):
        i = np.where(np.array(['#'.join(str(i) for i in z)
                               for z in zip(mut_array[:, 1], mut_array[:, 2])]) == key1)[0][0]
        ref = mut_array[i, 9]
        alt = mut_array[i, 10]
        pure_tags_dict[key1] = {}
        used_tags = []
        for key2, value2 in sorted(value1.items()):
            for key3, value3 in value2.items():
                pure_tag = key2[:-5]
                if key3 == alt:
                    if pure_tag in pure_tags_dict[key1]:
                        pure_tags_dict[key1][pure_tag] += 1
                    else:
                        pure_tags_dict[key1][pure_tag] = 1

    # 6. create pure_tags_dict_short with thresh
    if thresh > 0:
        pure_tags_dict_short = {}
        for key, value in sorted(pure_tags_dict.items()):
            if len(value) < thresh:
                    pure_tags_dict_short[key] = value
    else:
        pure_tags_dict_short = pure_tags_dict

    whole_array = []
    for k in pure_tags_dict.values():
        if len(k) != 0:
            keys = k.keys()
            if len(keys) > 1:
                for k1 in keys:
                    whole_array.append(k1)
            else:
                whole_array.append(keys[0])

    # 7. output summary with threshold
    workbook = xlsxwriter.Workbook(outfile)
    ws1 = workbook.add_worksheet("Results")
    ws2 = workbook.add_worksheet("Allele frequencies")
    ws3 = workbook.add_worksheet("Tiers")

    format1 = workbook.add_format({'bg_color': '#BCF5A9'})  # green
    format2 = workbook.add_format({'bg_color': '#FFC7CE'})  # red
    format3 = workbook.add_format({'bg_color': '#FACC2E'})  # yellow

    header_line = ('tag', 'mate', 'chr', 'pos', 'read pos.ab', 'read pos.ba', 'read median length.ab',
                   'read median length.ba', 'DCS median length', 'ref', 'alt',
                   'FS.ab', 'FS.ba' , 'FSqc.ab', 'FSqc.ba', 'ref.ab', 'ref.ba', 'alt.ab', 'alt.ba' ,
                   'rel. ref.ab', 'rel. ref.ba', 'rel. alt.ab', 'rel. alt.ba',
                   'na.ab', 'na.ba' , 'lowq.ab', 'lowq.ba',
                   'SSCS alt.ab', 'SSCS alt.ba', 'SSCS ref.ab', 'SSCS ref.ba',
                   'other mut', 'tiers', 'chimeric tag', 'variant ID')
    ws1.write_row(0, 0, header_line)

    counter_tier11 = 0
    counter_tier12 = 0
    counter_tier21 = 0
    counter_tier22 = 0
    counter_tier23 = 0
    counter_tier24 = 0
    counter_tier31 = 0
    counter_tier4 = 0

    row = 1
    tier_dict = {}
    for key1, value1 in sorted(mut_dict.items()):
        counts_mut = 0
        if key1 in pure_tags_dict_short.keys():
            i = np.where(np.array(['#'.join(str(i) for i in z)
                                   for z in zip(mut_array[:, 1], mut_array[:, 2])]) == key1)[0][0]
            ref = mut_array[i, 9]
            alt = mut_array[i, 10]
            dcs_median = cvrg_dict[key1][2]
            
            tier_dict[key1] = {}
            values_tier_dict = [("tier 1.1", 0), ("tier 1.2", 0), ("tier 2.1", 0), ("tier 2.2", 0),
                        ("tier 2.3", 0), ("tier 2.4", 0), ("tier 3.1", 0), ("tier 4", 0)]
            for k, v in values_tier_dict:
                tier_dict[key1][k] = v
            
            used_keys = []
            if 'ab' in mut_pos_dict[key1].keys():
                sscs_mut_ab = mut_pos_dict[key1]['ab']
            else:
                sscs_mut_ab = 0
            if 'ba' in mut_pos_dict[key1].keys():
                sscs_mut_ba = mut_pos_dict[key1]['ba']
            else:
                sscs_mut_ba = 0
            if 'ab' in ref_pos_dict[key1].keys():
                sscs_ref_ab = ref_pos_dict[key1]['ab']
            else:
                sscs_ref_ab = 0
            if 'ba' in ref_pos_dict[key1].keys():
                sscs_ref_ba = ref_pos_dict[key1]['ba']
            else:
                sscs_ref_ba = 0
            for key2, value2 in sorted(value1.items()):

                add_mut14 = ""
                add_mut23 = ""
                if (key2[:-5] in pure_tags_dict_short[key1].keys()) and (key2[:-5] not in used_keys):
                    if key2[:-5] + '.ab.1' in mut_dict[key1].keys():
                        total1 = sum(mut_dict[key1][key2[:-5] + '.ab.1'].values())
                        if 'na' in mut_dict[key1][key2[:-5] + '.ab.1'].keys():
                            na1 = mut_dict[key1][key2[:-5] + '.ab.1']['na']
                            na1f = na1/total1
                        else:
                            na1 = na1f = 0
                        if 'lowQ' in mut_dict[key1][key2[:-5] + '.ab.1'].keys():
                            lowq1 = mut_dict[key1][key2[:-5] + '.ab.1']['lowQ']
                            lowq1f = lowq1/total1
                        else:
                            lowq1 = lowq1f = 0
                        if ref in mut_dict[key1][key2[:-5] + '.ab.1'].keys():
                            ref1 = mut_dict[key1][key2[:-5] + '.ab.1'][ref]
                            ref1f = ref1/(total1-na1-lowq1)
                        else:
                            ref1 = ref1f = 0
                        if alt in mut_dict[key1][key2[:-5] + '.ab.1'].keys():
                            alt1 = mut_dict[key1][key2[:-5] + '.ab.1'][alt]
                            alt1f = alt1 / (total1 - na1 - lowq1)
                        else:
                            alt1 = alt1f = 0
                        total1new = total1 - na1 - lowq1
                        if (key2[:-5] + '.ab.1') in mut_read_dict.keys():
                            k1 = mut_read_dict[(key2[:-5] + '.ab.1')].keys()
                            add_mut1 = len(k1)
                            if add_mut1 > 1:
                                for k, v in mut_read_dict[(key2[:-5] + '.ab.1')].items():
                                    if k != key1:
                                        if len(add_mut14) == 0:
                                            add_mut14 = str(k) + "_" + v
                                        else:
                                            add_mut14 = add_mut14 + ", " + str(k) + "_" + v
                        else:
                            k1 = []
                    else:
                        total1 = total1new = na1 = lowq1 = 0
                        ref1 = alt1 = ref1f = alt1f = 0

                    if key2[:-5] + '.ab.2' in mut_dict[key1].keys():
                        total2 = sum(mut_dict[key1][key2[:-5] + '.ab.2'].values())
                        if 'na' in mut_dict[key1][key2[:-5] + '.ab.2'].keys():
                            na2 = mut_dict[key1][key2[:-5] + '.ab.2']['na']
                            na2f = na2 / total2
                        else:
                            na2 = na2f = 0
                        if 'lowQ' in mut_dict[key1][key2[:-5] + '.ab.2'].keys():
                            lowq2 = mut_dict[key1][key2[:-5] + '.ab.2']['lowQ']
                            lowq2f = lowq2 / total2
                        else:
                            lowq2 = lowq2f = 0
                        if ref in mut_dict[key1][key2[:-5] + '.ab.2'].keys():
                            ref2 = mut_dict[key1][key2[:-5] + '.ab.2'][ref]
                            ref2f = ref2 / (total2 - na2 - lowq2)
                        else:
                            ref2 = ref2f = 0
                        if alt in mut_dict[key1][key2[:-5] + '.ab.2'].keys():
                            alt2 = mut_dict[key1][key2[:-5] + '.ab.2'][alt]
                            alt2f = alt2 / (total2 - na2 - lowq2)
                        else:
                            alt2 = alt2f = 0
                        total2new = total2 - na2 - lowq2
                        if (key2[:-5] + '.ab.2') in mut_read_dict.keys():
                            k2 = mut_read_dict[(key2[:-5] + '.ab.2')].keys()
                            add_mut2 = len(k2)
                            if add_mut2 > 1:
                                for k, v in mut_read_dict[(key2[:-5] + '.ab.2')].items():
                                    if k != key1:
                                        if len(add_mut23) == 0:
                                            add_mut23 = str(k) + "_" + v
                                        else:
                                            add_mut23 = add_mut23 + ", " + str(k) + "_" + v
                        else:
                            k2 = []
                    else:
                        total2 = total2new = na2 = lowq2 = 0
                        ref2 = alt2 = ref2f = alt2f = 0

                    if key2[:-5] + '.ba.1' in mut_dict[key1].keys():
                        total3 = sum(mut_dict[key1][key2[:-5] + '.ba.1'].values())
                        if 'na' in mut_dict[key1][key2[:-5] + '.ba.1'].keys():
                            na3 = mut_dict[key1][key2[:-5] + '.ba.1']['na']
                            na3f = na3 / total3
                        else:
                            na3 = na3f = 0
                        if 'lowQ' in mut_dict[key1][key2[:-5] + '.ba.1'].keys():
                            lowq3 = mut_dict[key1][key2[:-5] + '.ba.1']['lowQ']
                            lowq3f = lowq3 / total3
                        else:
                            lowq3 = lowq3f = 0
                        if ref in mut_dict[key1][key2[:-5] + '.ba.1'].keys():
                            ref3 = mut_dict[key1][key2[:-5] + '.ba.1'][ref]
                            ref3f = ref3 / (total3 - na3 - lowq3)
                        else:
                            ref3 = ref3f = 0
                        if alt in mut_dict[key1][key2[:-5] + '.ba.1'].keys():
                            alt3 = mut_dict[key1][key2[:-5] + '.ba.1'][alt]
                            alt3f = alt3 / (total3 - na3 - lowq3)
                        else:
                            alt3 = alt3f = 0
                        total3new = total3 - na3 - lowq3
                        if (key2[:-5] + '.ba.1') in mut_read_dict.keys():
                            add_mut3 = len(mut_read_dict[(key2[:-5] + '.ba.1')].keys())
                            if add_mut3 > 1:
                                for k, v in mut_read_dict[(key2[:-5] + '.ba.1')].items():
                                    if k != key1 and k not in k2:
                                        if len(add_mut23) == 0:
                                            add_mut23 = str(k) + "_" + v
                                        else:
                                            add_mut23 = add_mut23 + ", " + str(k) + "_" + v
                    else:
                        total3 = total3new = na3 = lowq3 = 0
                        ref3 = alt3 = ref3f = alt3f = 0

                    if key2[:-5] + '.ba.2' in mut_dict[key1].keys():
                        total4 = sum(mut_dict[key1][key2[:-5] + '.ba.2'].values())
                        if 'na' in mut_dict[key1][key2[:-5] + '.ba.2'].keys():
                            na4 = mut_dict[key1][key2[:-5] + '.ba.2']['na']
                            na4f = na4 / total4
                        else:
                            na4 = na4f = 0
                        if 'lowQ' in mut_dict[key1][key2[:-5] + '.ba.2'].keys():
                            lowq4 = mut_dict[key1][key2[:-5] + '.ba.2']['lowQ']
                            lowq4f = lowq4 / total4
                        else:
                            lowq4 = lowq4f = 0
                        if ref in mut_dict[key1][key2[:-5] + '.ba.2'].keys():
                            ref4 = mut_dict[key1][key2[:-5] + '.ba.2'][ref]
                            ref4f = ref4 / (total4 - na4 - lowq4)
                        else:
                            ref4 = ref4f = 0
                        if alt in mut_dict[key1][key2[:-5] + '.ba.2'].keys():
                            alt4 = mut_dict[key1][key2[:-5] + '.ba.2'][alt]
                            alt4f = alt4 / (total4 - na4 - lowq4)
                        else:
                            alt4 = alt4f = 0
                        total4new = total4 - na4 - lowq4
                        if (key2[:-5] + '.ba.2') in mut_read_dict.keys():
                            add_mut4 = len(mut_read_dict[(key2[:-5] + '.ba.2')].keys())
                            if add_mut4 > 1:
                                for k, v in mut_read_dict[(key2[:-5] + '.ba.2')].items():
                                    if k != key1 and k not in k1:
                                        if len(add_mut14) == 0:
                                            add_mut14 = str(k) + "_" + v
                                        else:
                                            add_mut14 = add_mut14 + ", " + str(k) + "_" + v
                    else:
                        total4 = total4new = na4 = lowq4 = 0
                        ref4 = alt4 = ref4f = alt4f = 0

                    read_pos1 = read_pos2 = read_pos3 = read_pos4 = "None"
                    read_len_median1 = read_len_median2 = read_len_median3 = read_len_median4 = "None"
                    
                    if key2[:-5] + '.ab.1' in mut_read_pos_dict[key1].keys():
                        read_pos1 = np.median(mut_read_pos_dict[key1][key2[:-5] + '.ab.1'])
                        read_len_median1 = np.median(reads_dict[key1][key2[:-5] + '.ab.1'])
                    if key2[:-5] + '.ab.2' in mut_read_pos_dict[key1].keys():
                        read_pos2 = np.median(mut_read_pos_dict[key1][key2[:-5] + '.ab.2'])
                        read_len_median2 = np.median(reads_dict[key1][key2[:-5] + '.ab.2'])
                    if key2[:-5] + '.ba.1' in mut_read_pos_dict[key1].keys():
                        read_pos3 = np.median(mut_read_pos_dict[key1][key2[:-5] + '.ba.1'])
                        read_len_median3 = np.median(reads_dict[key1][key2[:-5] + '.ba.1'])
                    if key2[:-5] + '.ba.2' in mut_read_pos_dict[key1].keys():
                        read_pos4 = np.median(mut_read_pos_dict[key1][key2[:-5] + '.ba.2'])
                        read_len_median4 = np.median(reads_dict[key1][key2[:-5] + '.ba.2'])

                    used_keys.append(key2[:-5])
                    counts_mut += 1
                    if (alt1f + alt2f + alt3f + alt4f) > 0.5:
                        if total1new == 0:
                            ref1f = alt1f = "None"
                            alt1ff = -1
                        else:
                            alt1ff = alt1f
                        if total2new == 0:
                            ref2f = alt2f = "None"
                            alt2ff = -1
                        else:
                            alt2ff = alt2f
                        if total3new == 0:
                            ref3f = alt3f = "None"
                            alt3ff = -1
                        else:
                            alt3ff = alt3f
                        if total4new == 0:
                            ref4f = alt4f = "None"
                            alt4ff = -1
                        else:
                            alt4ff = alt4f
                        chrom, pos = re.split('\#', key1)
                        # assign tiers
                        if ((all(int(i) >= 3 for i in [total1new, total4new]) &
                             all(float(i) >= 0.75 for i in [alt1ff, alt4ff])) |
                                (all(int(i) >= 3 for i in [total2new, total3new]) &
                                 all(float(i) >= 0.75 for i in [alt2ff, alt3ff]))):
                            tier = "1.1"
                            counter_tier11 += 1
                            tier_dict[key1]["tier 1.1"] += 1

                        elif (all(int(i) >= 1 for i in [total1new, total2new, total3new, total4new]) &
                              any(int(i) >= 3 for i in [total1new, total4new]) &
                              any(int(i) >= 3 for i in [total2new, total3new]) &
                              all(float(i) >= 0.75 for i in [alt1ff, alt2ff, alt3ff, alt4ff])):
                            tier = "1.2"
                            counter_tier12 += 1
                            tier_dict[key1]["tier 1.2"] += 1

                        elif ((all(int(i) >= 1 for i in [total1new, total4new]) &
                               any(int(i) >= 3 for i in [total1new, total4new]) &
                               all(float(i) >= 0.75 for i in [alt1ff, alt4ff])) |
                              (all(int(i) >= 1 for i in [total2new, total3new]) &
                               any(int(i) >= 3 for i in [total2new, total3new]) &
                               all(float(i) >= 0.75 for i in [alt2ff, alt3ff]))):
                            tier = "2.1"
                            counter_tier21 += 1
                            tier_dict[key1]["tier 2.1"] += 1

                        elif (all(int(i) >= 1 for i in [total1new, total2new, total3new, total4new]) &
                              all(float(i) >= 0.75 for i in [alt1ff, alt2ff, alt3ff, alt4ff])):
                            tier = "2.2"
                            counter_tier22 += 1
                            tier_dict[key1]["tier 2.2"] += 1

                        elif ((all(int(i) >= 1 for i in [total1new, total4new]) &
                               any(int(i) >= 3 for i in [total2new, total3new]) &
                               all(float(i) >= 0.75 for i in [alt1ff, alt4ff]) &
                               any(float(i) >= 0.75 for i in [alt2ff, alt3ff])) |
                              (all(int(i) >= 1 for i in [total2new, total3new]) &
                               any(int(i) >= 3 for i in [total1new, total4new]) &
                               all(float(i) >= 0.75 for i in [alt2ff, alt3ff]) &
                               any(float(i) >= 0.75 for i in [alt1ff, alt4ff]))):
                            tier = "2.3"
                            counter_tier23 += 1
                            tier_dict[key1]["tier 2.3"] += 1

                        elif ((len(pure_tags_dict_short[key1]) > 1) &
                              (all(float(i) >= 0.5 for i in [alt1ff, alt4ff]) |
                               all(float(i) >= 0.5 for i in [alt2ff, alt3ff]))):
                            tier = "2.4"
                            counter_tier24 += 1
                            tier_dict[key1]["tier 2.4"] += 1

                        elif ((all(int(i) >= 1 for i in [total1new, total4new]) &
                               all(float(i) >= 0.5 for i in [alt1ff, alt4ff])) |
                              (all(int(i) >= 1 for i in [total2new, total3new]) &
                               all(float(i) >= 0.5 for i in [alt2ff, alt3ff]))):
                            tier = "3.1"
                            counter_tier31 += 1
                            tier_dict[key1]["tier 3.1"] += 1
                        else:
                            tier = "4"
                            counter_tier4 += 1
                            tier_dict[key1]["tier 4"] += 1

                        var_id = '-'.join([chrom, pos, ref, alt])
                        sample_tag = key2[:-5]
                        array2 = np.unique(whole_array)  # remove duplicate sequences to decrease running time
                        # exclude identical tag from array2, to prevent comparison to itself
                        same_tag = np.where(array2 == sample_tag)
                        index_array2 = np.arange(0, len(array2), 1)
                        index_withoutSame = np.delete(index_array2, same_tag)  # delete identical tag from the data
                        array2 = array2[index_withoutSame]
            if len(array2) != 0:  # only perform chimera analysis if there is more than 1 variant
                array1_half = sample_tag[0:int(len(sample_tag) / 2)]  # mate1 part1
                array1_half2 = sample_tag[int(len(sample_tag) / 2):int(len(sample_tag))]  # mate1 part 2
                array2_half = np.array([i[0:int(len(i) / 2)] for i in array2])  # mate2 part1
                array2_half2 = np.array([i[int(len(i) / 2):int(len(i))] for i in array2])  # mate2 part2

                min_tags_list_zeros = []
                chimera_tags = []
                for mate_b in [False, True]:
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
                    # calculate HD of "a" in the tag to all "a's" or "b" in the tag to all "b's"
                    dist = np.array([sum(itertools.imap(operator.ne, half1_mate1, c)) for c in half1_mate2])
                    min_index = np.where(dist == dist.min())  # get index of min HD
                    # get all "b's" of the tag or all "a's" of the tag with minimum HD
                    min_tag_half2 = half2_mate2[min_index]
                    min_tag_array2 = array2[min_index]  # get whole tag with min HD
                    min_value = dist.min()
                    # calculate HD of "b" to all "b's" or "a" to all "a's"
                    dist_second_half = np.array([sum(itertools.imap(operator.ne, half2_mate1, e))
                                                 for e in min_tag_half2])

                    dist2 = dist_second_half.max()
                    max_index = np.where(dist_second_half == dist_second_half.max())[0]  # get index of max HD
                    max_tag = min_tag_array2[max_index]

                    # tags which have identical parts:
                    if min_value == 0 or dist2 == 0:
                        min_tags_list_zeros.append(tag)
                        chimera_tags.append(max_tag)
                        chimeric = True
                    else:
                        chimeric = False

                    if mate_b is False:
                        text = "pos {}: sample tag: {}; HD a = {}; HD b' = {}; similar tag(s): {}; chimeric = {}". \
                                format(pos, sample_tag, min_value, dist2, list(max_tag), chimeric)
                    else:
                        text = "pos {}: sample tag: {}; HD a' = {}; HD b = {}; similar tag(s): {}; chimeric = {}". \
                                format(pos, sample_tag, dist2, min_value, list(max_tag), chimeric)
                    i += 1
                chimera_tags = [x for x in chimera_tags if x != []]
                chimera_tags_new = []
                for i in chimera_tags:
                    if len(i) > 1:
                        for t in i:
                            chimera_tags_new.append(t)
                    else:
                        chimera_tags_new.extend(i)
                chimera_tags_new = np.asarray(chimera_tags_new)
                chimera = ", ".join(chimera_tags_new)
            else:
                chimera = ""

            l = (key2[:-5], 'ab1.ba2', chrom, pos, read_pos1, read_pos4, read_len_median1,
                 read_len_median4, dcs_median, ref, alt,
                 total1, total4, total1new, total4new,
                 ref1, ref4, alt1, alt4, ref1f, ref4f, alt1f, alt4f,
                 na1, na4, lowq1, lowq4,
                 sscs_mut_ab, sscs_mut_ba, sscs_ref_ab, sscs_ref_ba,
                 add_mut14, tier, chimera, var_id)
            ws1.write_row(row, 0, l)
            l = (key2[:-5], 'ab2.ba1', chrom, pos, read_pos2, read_pos3, read_len_median2,
                 read_len_median3, dcs_median, ref, alt,
                 total2, total3, total2new, total3new,
                 ref2, ref3, alt2, alt3, ref2f, ref3f, alt2f, alt3f,
                 na2, na3, lowq2, lowq3,
                 sscs_mut_ab, sscs_mut_ba, sscs_ref_ab, sscs_ref_ba,
                 add_mut23, "", chimera, "")
            ws1.write_row(row + 1, 0, l)

            ws1.conditional_format('N{}:O{}'.format(row + 1, row + 2), 
                                   {'type': 'formula',
                                    'criteria': '=OR($AG${}="1.1", $AG${}="1.2")'.format(row + 1, row + 1),
                                    'format': format1, 
                                    'multi_range': 
                                        'N{}:O{} V{}:W{} AG{}'.format(row + 1, row + 2, row + 1, 
                                                                      row + 2, row + 1, row + 2)})

            ws1.conditional_format('N{}:O{}'.format(row + 1, row + 2),
                                   {'type': 'formula',
                                    'criteria': '=OR($AG${}="2.1",$AG${}="2.2", $AG${}="2.3", $AG${}="2.4")'.format(row + 1, row + 1, row + 1, row + 1),'format': format3, 'multi_range': 'N{}:O{} V{}:W{} AG{}'.format(row + 1, row + 2, row + 1, row + 2, row + 1, row + 2)})
            ws1.conditional_format('N{}:O{}'.format(row + 1, row + 2), 
                                   {'type': 'formula',
                                    'criteria': '=$AG${}>="3"'.format(row + 1),
                                    'format': format2, 
                                    'multi_range': 'N{}:O{} V{}:W{} AG{}'.format(row + 1, row + 2, row + 1, 
                                                                                 row + 2, row + 1, row + 2)})

            row += 3
   
    # sheet 2 
    header_line2 = ('variant ID', 'cvrg (Du Novo)', 'tier 1.1', 'AF 1.1', 'tier 1.2', 'AF 1.2',
                    'tier 2.1', 'AF 2.1', 'tier 2.2', 'AF 2.2', 'tier 2.3', 'AF 2.3', 'tier 2.4', 'AF 2.4',
                    'tier 3.1', 'AF 3.1', 'tier 4', 'AF 4', 'AF 1.1-1.2', 'AF 1.1-2.1', 'AF 1.1-2.2', 
                    'AF 1.1-2.3', 'AF 1.1-2.4', 'AF 1.1-3.1', 'AF 1.1-4', 'cvrg (VA)', 'tier count', 
                    'AF (VA)', 'alt. allele (Du Novo)', 'AF (Du Novo)', 'ref. allele (Du Novo)')
    ws2.write_row(0, 0, header_line2)
    row = 0

    for key1, value1 in sorted(tier_dict.items()):
        if key1 in pure_tags_dict_short.keys():
            i = np.where(np.array(['#'.join(str(i) for i in z) 
                                   for z in zip(mut_array[:, 1], mut_array[:, 2])]) == key1)[0][0]
            ref = mut_array[i, 9]
            alt = mut_array[i, 10]
            cvrg = mut_array[i, 7]
            chrom, pos = re.split('\#', key1)
            ref_count = cvrg_dict[key1][0]
            alt_count = cvrg_dict[key1][1]

            var_id = '-'.join([chrom, pos, ref, alt])
            l = [var_id, ref_count+alt_count]
            used_tiers = []
            cum_af = []
            for key2, value2 in sorted(value1.items()):
                l.append(value2)
                if (value2 + ref_count) != 0:
                    af = value2 / (value2 + ref_count)
                else:
                    af = 0
                l.append(af)

                # calculate cummulative AF
                used_tiers.append(value2)
                if len(used_tiers) > 1:
                    if (sum(used_tiers) + ref_count) != 0:
                        cum = sum(used_tiers) / (sum(used_tiers) + ref_count)
                    else:
                        cum = 0
                    cum_af.append(cum)
            l.extend(cum_af)
            l.extend([ref_count + sum(used_tiers), sum(used_tiers), sum(used_tiers) / (ref_count + sum(used_tiers)), 
                      alt_count, alt_count / (ref_count + alt_count), ref_count])
            
            l = tuple(l)
            ws2.write_row(row + 1, 0, l)
            ws2.conditional_format('C{}:F{}'.format(row + 2, row + 2),  {'type': 'formula',
                                            'criteria': '=$C$1="tier 1.1"',
                                            'format': format1,
                                            'multi_range': 'C{}:F{} C1:F1'.format(row + 2, row + 2)})
            ws2.conditional_format('G{}:N{}'.format(row + 2, row + 2),  {'type': 'formula',
                                            'criteria': '=$G$1="tier 2.1"',
                                            'format': format3,
                                            'multi_range': 'G{}:N{} G1:N1'.format(row + 2, row + 2)})
            ws2.conditional_format('O{}:R{}'.format(row + 2, row + 2),  {'type': 'formula',
                                            'criteria': '=$O$1="tier 3.1"',
                                            'format': format2,
                                            'multi_range': 'O{}:R{} O1:R1'.format(row + 2, row + 2)})
            row += 1
            
    # sheet 3
    sheet3 = [("tier 1.1", counter_tier11), ("tier 1.2", counter_tier12), ("tier 2.1", counter_tier21),
              ("tier 2.2", counter_tier22), ("tier 2.3", counter_tier23), ("tier 2.4", counter_tier24),
              ("tier 3.1", counter_tier31), ("tier 4", counter_tier4)]

    header = ("tier", "count")
    ws3.write_row(0, 0, header)
    
    for i in range(len(sheet3)):
        ws3.write_row(i+1, 0, sheet3[i])
        ws3.conditional_format('A{}:B{}'.format(i+2, i+2), 
                               {'type': 'formula',
                                'criteria': '=OR($A${}="tier 1.1", $A${}="tier 1.2")'.format(i+2, i+2),
                                'format': format1})
        ws3.conditional_format('A{}:B{}'.format(i+2, i+2), 
                               {'type': 'formula',
                                'criteria': '=OR($A${}="tier 2.1", $A${}="tier 2.2", $A${}="tier 2.3", $A${}="tier 2.4")'.format(i+2, i+2, i+2, i+2),
                                'format': format3})
        ws3.conditional_format('A{}:B{}'.format(i+2, i+2), 
                               {'type': 'formula',
                                'criteria': '=OR($A${}="tier 3.1", $A${}="tier 4")'.format(i+2, i+2),
                                'format': format2})

    description_tiers = [("Tier 1.1", "both ab and ba SSCS present (>75% of the sites with alternative base) and minimal FS>=3 for both SSCS in at least one mate"),
            ("", ""),
            ("Tier 1.2", "both ab and ba SSCS present (>75% of the sites with alt. base) and mate pair validation (min. FS=1) and minimal FS>=3 for at least one of the SSCS"),
            ("Tier 2.1", "both ab and ba SSCS present (>75% of the sites with alt. base) and minimal FS>=3 for at least one of the SSCS in at least one mate"),
            ("Tier 2.2", "both ab and ba SSCS present (>75% of the sites with alt. base) and mate pair validation (min. FS=1)"),
            ("Tier 2.3", "both ab and ba SSCS present (>75% of the sites with alt. base) and minimal FS=1 for both SSCS in one mate and minimal FS>=3 for at least one of the SSCS in the other mate"),
            ("Tier 2.4", "both ab and ba SSCS present (>50% of the sites with alt. base) and recurring mutation on this position"),
            ("Tier 3.1", "both ab and ba SSCS present (>50% of the sites with alt. base) and minimal FS>=1 for both SSCS in at least one mate"),
            ("Tier 4", "remaining variants")]
    examples_tiers = [[("AAAAAGATGCCGACTACCTT", "ab1.ba2", "Chr5:5-20000", "11068", "254", "228", "287", "288", "289", 
                        "C", "G", "3", "6", "3", "6", "0", "0", "3", "6", "0", "0", "1", "1", "0", "0", "0", "0", 
                        "4081", "4098", "5", "10", "", "1.1", "", "Chr5:5-20000-11068-C-G"), 
                       ("AAAAAGATGCCGACTACCTT", "ab2.ba1", "Chr5:5-20000", "11068", "None", "None", "None", "None", 
                        "289", "C", "G", "0", "0", "0", "0", "0", "0", "3", "6", "None", "None", "None", "None", 
                        "0", "0", "0", "0", "4081", "4098", "5", "10", "", "", "", "")],
                        [("AAAAATGCGTAGAAATATGC", "ab1.ba2", "Chr5:5-20000", "11068", "254", "228", "287", "288", "289", 
                          "C", "G", "33", "43", "33", "43", "0", "0", "33", "43", "0", "0", "1", "1", "0", "0", "0", 
                          "0", "4081", "4098", "5", "10", "", "1.1", "", "Chr5:5-20000-11068-C-G"), 
                         ("AAAAATGCGTAGAAATATGC", "ab2.ba1", "Chr5:5-20000", "11068", "268", "268", "270", "288", "289", 
                          "C", "G", "11", "34", "10", "27", "0", "0", "10", "27", "0", "0", "1", "1", "0", "0", "1", 
                          "7", "4081", "4098", "5", "10", "", "", "", "")],
                        [("CTATGACCCGTGAGCCCATG", "ab1.ba2", "Chr5:5-20000", "10776", "132", "132", "287", "288", "290", 
                          "G", "T", "4", "1", "4", "1", "0", "0", "4", "1", "0", "0", "1", "1", "0", "0", "0", "0", "1", 
                          "6", "47170", "41149", "", "1.2", "", "Chr5:5-20000-10776-G-T"), 
                         ("CTATGACCCGTGAGCCCATG", "ab2.ba1", "Chr5:5-20000", "10776", "77", "132", "233", "200", "290", 
                          "G", "T", "4", "1", "4", "1", "0", "0", "4", "1", "0", "0", "1", "1", "0", "0", "0", "0", "1", 
                          "6", "47170", "41149", "", "", "", "")],
                        [("AAAAAAACATCATACACCCA", "ab1.ba2", "Chr5:5-20000", "11068", "246", "244", "287", "288", "289", 
                          "C", "G", "2", "8", "2", "8", "0", "0", "2", "8", "0", "0", "1", "1", "0", "0", "0", "0", 
                          "4081", "4098", "5", "10", "", "2.1", "", "Chr5:5-20000-11068-C-G"), 
                         ("AAAAAAACATCATACACCCA", "ab2.ba1", "Chr5:5-20000", "11068", "None", "None", "None", "None", 
                          "289", "C", "G", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "None", "None", "0", "0", 
                          "0", "0", "4081", "4098", "5", "10", "", "", "", "")],
                        [("ATCAGCCATGGCTATTATTG", "ab1.ba2", "chrom5:5-20000", "11068", "72", "72", "217", "288", "289",
                          "C", "G", "1", "1", "1", "1", "0", "0", "1", "1", "0", "0", "1", "1", "0", "0", "0", "0", 
                          "4081", "4098", "5", "10", "", "2.2", "", "Chr5:5-20000-11068-C-G"),
                         ("ATCAGCCATGGCTATTATTG", "ab2.ba1", "Chr5:5-20000", "11068", "153", "164", "217", "260", "289", 
                          "C", "G", "1", "1", "1", "1", "0", "0", "1", "1", "0", "0", "1", "1", "0", "0", "0", "0", 
                          "4081", "4098", "5", "10", "", "", "", "")],
                        [("ATCAATATGGCCTCGCCACG", "ab1.ba2", "Chr5:5-20000", "11068", "None", "None", "None", "None", 
                          "289", "C", "G", "0", "5", "0", "5", "0", "0", "0", "5", "None", "None", "None", "1", "0", 
                          "0", "0", "0", "4081", "4098", "5", "10", "", "2.3", "", "Chr5:5-20000-11068-C-G"),
                         ("ATCAATATGGCCTCGCCACG", "ab2.ba1", "Chr5:5-20000", "11068", "202", "255", "277", "290", "289", 
                          "C", "G", "1", "3", "1", "3", "0", "0", "1", "3", "0", "0", "1", "1", "0", "0", "1", "7", 
                          "4081", "4098", "5", "10", "", "", "", "")],
                        [("ATGCCTACCTCATTTGTCGT", "ab1.ba2", "Chr5:5-20000", "10776", "46", "15", "287", "288", "290", 
                          "G", "T", "3", "3", "3", "2", "3", "1", "0", "1", "1", "0.5", "0", "0.5", "0", "0", "0", "1", 
                          "3", "3", "47170", "41149", "", "2.4", "", "Chr5:5-20000-10776-G-T"),
                         ("ATGCCTACCTCATTTGTCGT", "ab2.ba1", "Chr5:5-20000", "10776", "None", "274", "None", 
                          "288", "290", "G", "T", "0", "3", "0", "2", "0", "1", "0", "1", "None", "0.5", "None", "0.5", 
                          "0", "0", "0", "1", "3", "3", "47170", "41149", "", "", "", "")],
                        [("ACAACATCACGTATTCAGGT", "ab1.ba2", "Chr5:5-20000", "11315", "197", "197", "240", "255", "271", 
                          "C", "T", "2", "3", "2", "3", "0", "1", "2", "2", "0", "0.333333333333333", "1", 
                          "0.666666666666667", "0", "0", "0", "0", "1", "1", "6584", "6482", "", "3.1", "", 
                          "Chr5:5-20000-11315-C-T"),
                         ("ACAACATCACGTATTCAGGT", "ab2.ba1", "Chr5:5-20000", "11315", "35", "35", "240", "258", "271", 
                          "C", "T", "2", "3", "2", "3", "0", "1", "2", "2", "0", "0.333333333333333", "1", 
                          "0.666666666666667", "0", "0", "0", "0", "1", "1", "6584", "6482", "", "", "", "")],
                        [("ATGTTGTGAATAACCCACAC", "ab1.ba2", "Chr5:5-20000", "13983", "209", "186", "255", "276", "269", 
                          "G", "C", "0", "6", "0", "6", "0", "0", "0", "6", "0", "0", "0", "1", "0", "0", "0", "0", "1", 
                          "1", "5348", "5350", "", "4", "", "Chr5:5-20000-13983-G-C"),
                         ("ATGTTGTGAATAACCCACAC", "ab2.ba1", "Chr5:5-20000", "13983", "None", "None", "None", "None", 
                          "269", "G", "C", "0", "0", "0", "0", "0", "0", "0", "0", "None", "None", "None", "None", "0", 
                          "0", "0", "0", "1", "1", "5348", "5350", "", "", "", "")]]

    ws3.write(10, 0, "Description of tiers with examples")
    ws3.write_row(11, 0, header_line)
    row = 0
    for i in range(len(description_tiers)):
        ws3.write_row(12+row+i+1, 0, description_tiers[i])
        ex = examples_tiers[i]
    for k in range(len(ex)):
        ws3.write_row(12+row+i+k+2, 0, ex[k])
        ws3.conditional_format('N{}:O{}'.format(12+row+i+k+2, 12+row+i+k+3), 
                              {'type': 'formula',
                              'criteria': '=OR($AG${}="1.1", $AG${}="1.2")'.format(12+row+i+k+2, 12+row+i+k+2),
                              'format': format1, 
                              'multi_range': 
                                   'N{}:O{} V{}:W{} AG{}'.format(12+row+i+k+2, 12+row+i+k+3, 12+row+i+k+2, 
                                                                 12+row+i+k+3, 12+row+i+k+2, 12+row+i+k+3)})

        ws3.conditional_format('N{}:O{}'.format(12+row+i+k+2, 12+row+i+k+3), 
                              {'type': 'formula',
                               'criteria': 
                                   '=OR($AG${}="2.1",$AG${}="2.2", $AG${}="2.3", $AG${}="2.4")'.format(12+row+i+k+2,
                                                                                                       12+row+i+k+2,
                                                                                                       12+row+i+k+2,
                                                                                                       12+row+i+k+2),
                               'format': format3, 
                               'multi_range': 
                                   'N{}:O{} V{}:W{} AG{}'.format(12+row+i+k+2, 12+row+i+k+3, 12+row+i+k+2, 
                                                                 12+row+i+k+3, 12+row+i+k+2, 12+row+i+k+3)})
        ws3.conditional_format('N{}:O{}'.format(12+row+i+k+2, 12+row+i+k+3), 
                               {'type': 'formula',
                                'criteria': '=$AG${}>="3"'.format(12+row+i+k+2),
                                'format': format2, 
                                'multi_range': 
                                    'N{}:O{} V{}:W{} AG{}'.format(12+row+i+k+2, 12+row+i+k+3, 12+row+i+k+2, 
                                                                  12+row+i+k+3, 12+row+i+k+2, 12+row+i+k+3)})
        row += 3
    workbook.close()


if __name__ == '__main__':
    sys.exit(read2mut(sys.argv))
