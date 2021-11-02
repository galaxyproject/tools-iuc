#!/usr/bin/env python

import argparse
import gzip
import os
from collections import OrderedDict

import yaml
from Bio.SeqIO.QualityIO import FastqGeneralIterator

OUTPUT_DBKEY_DIR = 'output_dbkey'
OUTPUT_METRICS_DIR = 'output_metrics'


def get_sample_name(file_path):
    base_file_name = os.path.basename(file_path)
    if base_file_name.find(".") > 0:
        # Eliminate the extension.
        return os.path.splitext(base_file_name)[0]
    return base_file_name


def get_dbkey(dnaprints_dict, key, s):
    # dnaprints_dict looks something like this:
    # {'brucella': {'NC_002945v4': ['11001110', '11011110', '11001100']}
    # {'bovis': {'NC_006895': ['11111110', '00010010', '01111011']}}
    d = dnaprints_dict.get(key, {})
    for data_table_value, v_list in d.items():
        if s in v_list:
            return data_table_value
    return ""


def get_dnaprints_dict(dnaprint_fields):
    # A dndprint_fields entry looks something liek this.
    # [['AF2122', '/galaxy/tool-data/vsnp/AF2122/dnaprints/NC_002945v4.yml']]
    dnaprints_dict = {}
    for item in dnaprint_fields:
        # Here item is a 2-element list of data
        # table components, # value and path.
        value = item[0]
        path = item[1].strip()
        with open(path, "rt") as fh:
            # The format of all dnaprints yaml
            # files is something like this:
            # brucella:
            #  - 0111111111111111
            print_dict = yaml.load(fh, Loader=yaml.Loader)
        for print_dict_k, print_dict_v in print_dict.items():
            dnaprints_v_dict = dnaprints_dict.get(print_dict_k, {})
            if len(dnaprints_v_dict) > 0:
                # dnaprints_dict already contains k (e.g., 'brucella',
                # and dnaprints_v_dict will be a dictionary # that
                # looks something like this:
                # {'NC_002945v4': ['11001110', '11011110', '11001100']}
                value_list = dnaprints_v_dict.get(value, [])
                value_list = value_list + print_dict_v
                dnaprints_v_dict[value] = value_list
            else:
                # dnaprints_v_dict is an empty dictionary.
                dnaprints_v_dict[value] = print_dict_v
            dnaprints_dict[print_dict_k] = dnaprints_v_dict
    # dnaprints_dict looks something like this:
    # {'brucella': {'NC_002945v4': ['11001110', '11011110', '11001100']}
    # {'bovis': {'NC_006895': ['11111110', '00010010', '01111011']}}
    return dnaprints_dict


def get_group_and_dbkey(dnaprints_dict, brucella_string, brucella_sum, bovis_string, bovis_sum, para_string, para_sum):
    if brucella_sum > 3:
        group = "Brucella"
        dbkey = get_dbkey(dnaprints_dict, "brucella", brucella_string)
    elif bovis_sum > 3:
        group = "TB"
        dbkey = get_dbkey(dnaprints_dict, "bovis", bovis_string)
    elif para_sum >= 1:
        group = "paraTB"
        dbkey = get_dbkey(dnaprints_dict, "para", para_string)
    else:
        group = ""
        dbkey = ""
    return group, dbkey


def get_oligo_dict():
    oligo_dict = {}
    oligo_dict["01_ab1"] = "AATTGTCGGATAGCCTGGCGATAACGACGC"
    oligo_dict["02_ab3"] = "CACACGCGGGCCGGAACTGCCGCAAATGAC"
    oligo_dict["03_ab5"] = "GCTGAAGCGGCAGACCGGCAGAACGAATAT"
    oligo_dict["04_mel"] = "TGTCGCGCGTCAAGCGGCGTGAAATCTCTG"
    oligo_dict["05_suis1"] = "TGCGTTGCCGTGAAGCTTAATTCGGCTGAT"
    oligo_dict["06_suis2"] = "GGCAATCATGCGCAGGGCTTTGCATTCGTC"
    oligo_dict["07_suis3"] = "CAAGGCAGATGCACATAATCCGGCGACCCG"
    oligo_dict["08_ceti1"] = "GTGAATATAGGGTGAATTGATCTTCAGCCG"
    oligo_dict["09_ceti2"] = "TTACAAGCAGGCCTATGAGCGCGGCGTGAA"
    oligo_dict["10_canis4"] = "CTGCTACATAAAGCACCCGGCGACCGAGTT"
    oligo_dict["11_canis"] = "ATCGTTTTGCGGCATATCGCTGACCACAGC"
    oligo_dict["12_ovis"] = "CACTCAATCTTCTCTACGGGCGTGGTATCC"
    oligo_dict["13_ether2"] = "CGAAATCGTGGTGAAGGACGGGACCGAACC"
    oligo_dict["14_63B1"] = "CCTGTTTAAAAGAATCGTCGGAACCGCTCT"
    oligo_dict["15_16M0"] = "TCCCGCCGCCATGCCGCCGAAAGTCGCCGT"
    oligo_dict["16_mel1b"] = "TCTGTCCAAACCCCGTGACCGAACAATAGA"
    oligo_dict["17_tb157"] = "CTCTTCGTATACCGTTCCGTCGTCACCATGGTCCT"
    oligo_dict["18_tb7"] = "TCACGCAGCCAACGATATTCGTGTACCGCGACGGT"
    oligo_dict["19_tbbov"] = "CTGGGCGACCCGGCCGACCTGCACACCGCGCATCA"
    oligo_dict["20_tb5"] = "CCGTGGTGGCGTATCGGGCCCCTGGATCGCGCCCT"
    oligo_dict["21_tb2"] = "ATGTCTGCGTAAAGAAGTTCCATGTCCGGGAAGTA"
    oligo_dict["22_tb3"] = "GAAGACCTTGATGCCGATCTGGGTGTCGATCTTGA"
    oligo_dict["23_tb4"] = "CGGTGTTGAAGGGTCCCCCGTTCCAGAAGCCGGTG"
    oligo_dict["24_tb6"] = "ACGGTGATTCGGGTGGTCGACACCGATGGTTCAGA"
    oligo_dict["25_para"] = "CCTTTCTTGAAGGGTGTTCG"
    oligo_dict["26_para_sheep"] = "CGTGGTGGCGACGGCGGCGGGCCTGTCTAT"
    oligo_dict["27_para_cattle"] = "TCTCCTCGGTCGGTGATTCGGGGGCGCGGT"
    return oligo_dict


def get_seq_counts(value, fastq_list, gzipped):
    count = 0
    for fastq_file in fastq_list:
        if gzipped:
            with gzip.open(fastq_file, 'rt') as fh:
                for title, seq, qual in FastqGeneralIterator(fh):
                    count += seq.count(value)
        else:
            with open(fastq_file, 'r') as fh:
                for title, seq, qual in FastqGeneralIterator(fh):
                    count += seq.count(value)
    return(value, count)


def get_species_counts(fastq_list, gzipped):
    count_summary = {}
    oligo_dict = get_oligo_dict()
    for v1 in oligo_dict.values():
        returned_value, count = get_seq_counts(v1, fastq_list, gzipped)
        for key, v2 in oligo_dict.items():
            if returned_value == v2:
                count_summary.update({key: count})
    count_list = []
    for v in count_summary.values():
        count_list.append(v)
    brucella_sum = sum(count_list[:16])
    bovis_sum = sum(count_list[16:24])
    para_sum = sum(count_list[24:])
    return count_summary, count_list, brucella_sum, bovis_sum, para_sum


def get_species_strings(count_summary):
    binary_dictionary = {}
    for k, v in count_summary.items():
        if v > 1:
            binary_dictionary.update({k: 1})
        else:
            binary_dictionary.update({k: 0})
    binary_dictionary = OrderedDict(sorted(binary_dictionary.items()))
    binary_list = []
    for v in binary_dictionary.values():
        binary_list.append(v)
    brucella_binary = binary_list[:16]
    brucella_string = ''.join(str(e) for e in brucella_binary)
    bovis_binary = binary_list[16:24]
    bovis_string = ''.join(str(e) for e in bovis_binary)
    para_binary = binary_list[24:]
    para_string = ''.join(str(e) for e in para_binary)
    return brucella_string, bovis_string, para_string


def output_dbkey(file_name, dbkey, output_file):
    # Output the dbkey.
    with open(output_file, "w") as fh:
        fh.write("%s" % dbkey)


def output_files(fastq_file, count_list, group, dbkey, dbkey_file, metrics_file):
    base_file_name = get_sample_name(fastq_file)
    output_dbkey(base_file_name, dbkey, dbkey_file)
    output_metrics(base_file_name, count_list, group, dbkey, metrics_file)


def output_metrics(file_name, count_list, group, dbkey, output_file):
    # Output the metrics.
    with open(output_file, "w") as fh:
        fh.write("Sample: %s\n" % file_name)
        fh.write("Brucella counts: ")
        for i in count_list[:16]:
            fh.write("%d," % i)
        fh.write("\nTB counts: ")
        for i in count_list[16:24]:
            fh.write("%d," % i)
        fh.write("\nPara counts: ")
        for i in count_list[24:]:
            fh.write("%d," % i)
        fh.write("\nGroup: %s" % group)
        fh.write("\ndbkey: %s\n" % dbkey)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--dnaprint_fields', action='append', dest='dnaprint_fields', nargs=2, help="List of dnaprints data table value, name and path fields")
    parser.add_argument('--read1', action='store', dest='read1', help='Required: single read')
    parser.add_argument('--read2', action='store', dest='read2', required=False, default=None, help='Optional: paired read')
    parser.add_argument('--gzipped', action='store_true', dest='gzipped', help='Input files are gzipped')
    parser.add_argument('--output_dbkey', action='store', dest='output_dbkey', help='Output reference file')
    parser.add_argument('--output_metrics', action='store', dest='output_metrics', help='Output metrics file')

    args = parser.parse_args()

    fastq_list = [args.read1]
    if args.read2 is not None:
        fastq_list.append(args.read2)

    # The value of dnaprint_fields is a list of lists, where each list is
    # the [value, name, path] components of the vsnp_dnaprints data table.
    # The data_manager_vsnp_dnaprints tool assigns the dbkey column from the
    # all_fasta data table to the value column in the vsnp_dnaprints data
    # table to ensure a proper mapping for discovering the dbkey.
    dnaprints_dict = get_dnaprints_dict(args.dnaprint_fields)

    # Here fastq_list consists of either a single read
    # or a set of paired reads, producing single outputs.
    count_summary, count_list, brucella_sum, bovis_sum, para_sum = get_species_counts(fastq_list, args.gzipped)
    brucella_string, bovis_string, para_string = get_species_strings(count_summary)
    group, dbkey = get_group_and_dbkey(dnaprints_dict, brucella_string, brucella_sum, bovis_string, bovis_sum, para_string, para_sum)
    output_files(args.read1, count_list, group, dbkey, dbkey_file=args.output_dbkey, metrics_file=args.output_metrics)
