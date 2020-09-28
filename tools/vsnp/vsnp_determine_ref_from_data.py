#!/usr/bin/env python

import argparse
import gzip
import multiprocessing
import os
import queue
from collections import OrderedDict

import yaml
from Bio.SeqIO.QualityIO import FastqGeneralIterator

INPUT_READS_DIR = 'input_reads'
OUTPUT_DBKEY_DIR = 'output_dbkey'
OUTPUT_METRICS_DIR = 'output_metrics'


def get_base_file_name(file_path):
    base_file_name = os.path.basename(file_path)
    if base_file_name.find(".") > 0:
        # Eliminate the extension.
        return os.path.splitext(base_file_name)[0]
    elif base_file_name.find("_fq") > 0:
        # The "." character has likely
        # changed to an "_" character.
        return base_file_name.split("_fq")[0]
    elif base_file_name.find("_fastq") > 0:
        return base_file_name.split("_fastq")[0]
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


def get_group_and_dbkey_for_collection(task_queue, finished_queue, dnaprints_dict, timeout):
    while True:
        try:
            tup = task_queue.get(block=True, timeout=timeout)
        except queue.Empty:
            break
        fastq_file, count_list, brucella_string, brucella_sum, bovis_string, bovis_sum, para_string, para_sum = tup
        group, dbkey = get_group_and_dbkey(dnaprints_dict, brucella_string, brucella_sum, bovis_string, bovis_sum, para_string, para_sum)
        finished_queue.put((fastq_file, count_list, group, dbkey))
        task_queue.task_done()


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
        if gzipped == "true":
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


def get_species_counts_for_collection(task_queue, finished_queue, gzipped, timeout):
    while True:
        try:
            fastq_file = task_queue.get(block=True, timeout=timeout)
        except queue.Empty:
            break
        count_summary, count_list, brucella_sum, bovis_sum, para_sum = get_species_counts([fastq_file], gzipped)
        finished_queue.put((fastq_file, count_summary, count_list, brucella_sum, bovis_sum, para_sum))
        task_queue.task_done()


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


def get_species_strings_for_collection(task_queue, finished_queue, timeout):
    while True:
        try:
            tup = task_queue.get(block=True, timeout=timeout)
        except queue.Empty:
            break
        fastq_file, count_summary, count_list, brucella_sum, bovis_sum, para_sum = tup
        brucella_string, bovis_string, para_string = get_species_strings(count_summary)
        finished_queue.put((fastq_file, count_list, brucella_string, brucella_sum, bovis_string, bovis_sum, para_string, para_sum))
        task_queue.task_done()


def output_dbkey(file_name, dbkey, output_file=None):
    # Output the dbkey.
    if output_file is None:
        # We're producing a dataset collection.
        output_file = os.path.join(OUTPUT_DBKEY_DIR, "%s.txt" % file_name)
    with open(output_file, "w") as fh:
        fh.write("%s" % dbkey)


def output_files(fastq_file, count_list, group, dbkey, dbkey_file=None, metrics_file=None):
    base_file_name = get_base_file_name(fastq_file)
    output_dbkey(base_file_name, dbkey, dbkey_file)
    output_metrics(base_file_name, count_list, group, dbkey, metrics_file)


def output_files_for_collection(task_queue, timeout):
    while True:
        try:
            tup = task_queue.get(block=True, timeout=timeout)
        except queue.Empty:
            break
        fastq_file, count_list, group, dbkey = tup
        output_files(fastq_file, count_list, group, dbkey)
        task_queue.task_done()


def output_metrics(file_name, count_list, group, dbkey, output_file=None):
    # Output the metrics.
    if output_file is None:
        # We're producing a dataset collection.
        output_file = os.path.join(OUTPUT_METRICS_DIR, "%s.txt" % file_name)
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


def set_num_cpus(num_files, processes):
    num_cpus = int(multiprocessing.cpu_count())
    if num_files < num_cpus and num_files < processes:
        return num_files
    if num_cpus < processes:
        half_cpus = int(num_cpus / 2)
        if num_files < half_cpus:
            return num_files
        return half_cpus
    return processes


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--dnaprint_fields', action='append', dest='dnaprint_fields', nargs=2, required=False, default=None, help="List of dnaprints data table value, name and path fields")
    parser.add_argument('--read1', action='store', dest='read1', required=False, default=None, help='Required: single read')
    parser.add_argument('--read2', action='store', dest='read2', required=False, default=None, help='Optional: paired read')
    parser.add_argument('--gzipped', action='store', dest='gzipped', help='Input files are gzipped')
    parser.add_argument('--in_test_mode', action='store', dest='in_test_mode', required=False, default=None, help='Functional test mode flag')
    parser.add_argument('--output_dbkey', action='store', dest='output_dbkey', required=False, default=None, help='Output reference file')
    parser.add_argument('--output_metrics', action='store', dest='output_metrics', required=False, default=None, help='Output metrics file')
    parser.add_argument('--processes', action='store', dest='processes', type=int, help='User-selected number of processes to use for job splitting')

    args = parser.parse_args()

    collection = False
    fastq_list = []
    if args.read1 is not None:
        fastq_list.append(args.read1)
        if args.read2 is not None:
            fastq_list.append(args.read2)
    else:
        collection = True
        for file_name in sorted(os.listdir(INPUT_READS_DIR)):
            file_path = os.path.abspath(os.path.join(INPUT_READS_DIR, file_name))
            fastq_list.append(file_path)

    # The value of dnaprint_fields is a list of lists, where each list is
    # the [value, name, path] components of the vsnp_dnaprints data table.
    # The data_manager_vsnp_dnaprints tool assigns the dbkey column from the
    # all_fasta data table to the value column in the vsnp_dnaprints data
    # table to ensure a proper mapping for discovering the dbkey.
    if args.in_test_mode is None:
        dnaprints_dict = get_dnaprints_dict(args.dnaprint_fields)
    else:
        dnaprints_dict = {'bovis': {'AF2122': ['11001110', '11011110', '11001100']}}

    if collection:
        # Here fastq_list consists of any number of
        # reads, so each file will be processed and
        # dataset collections will be produced as outputs.
        multiprocessing.set_start_method('spawn')
        queue1 = multiprocessing.JoinableQueue()
        queue2 = multiprocessing.JoinableQueue()
        num_files = len(fastq_list)
        cpus = set_num_cpus(num_files, args.processes)
        # Set a timeout for get()s in the queue.
        timeout = 0.05

        for fastq_file in fastq_list:
            queue1.put(fastq_file)

        # Complete the get_species_counts task.
        processes = [multiprocessing.Process(target=get_species_counts_for_collection, args=(queue1, queue2, args.gzipped, timeout, )) for _ in range(cpus)]
        for p in processes:
            p.start()
        for p in processes:
            p.join()
        queue1.join()

        # Complete the get_species_strings task.
        processes = [multiprocessing.Process(target=get_species_strings_for_collection, args=(queue2, queue1, timeout, )) for _ in range(cpus)]
        for p in processes:
            p.start()
        for p in processes:
            p.join()
        queue2.join()

        # Complete the get_group_and_dbkey task.
        processes = [multiprocessing.Process(target=get_group_and_dbkey_for_collection, args=(queue1, queue2, dnaprints_dict, timeout, )) for _ in range(cpus)]
        for p in processes:
            p.start()
        for p in processes:
            p.join()
        queue1.join()

        # Complete the output_files task.
        processes = [multiprocessing.Process(target=output_files_for_collection, args=(queue2, timeout, )) for _ in range(cpus)]
        for p in processes:
            p.start()
        for p in processes:
            p.join()
        queue2.join()

        if queue1.empty() and queue2.empty():
            queue1.close()
            queue1.join_thread()
            queue2.close()
            queue2.join_thread()
    else:
        # Here fastq_list consists of either a single read
        # or a set of paired reads, producing single outputs.
        count_summary, count_list, brucella_sum, bovis_sum, para_sum = get_species_counts(fastq_list, args.gzipped)
        brucella_string, bovis_string, para_string = get_species_strings(count_summary)
        group, dbkey = get_group_and_dbkey(dnaprints_dict, brucella_string, brucella_sum, bovis_string, bovis_sum, para_string, para_sum)
        output_files(args.read1, count_list, group, dbkey, dbkey_file=args.output_dbkey, metrics_file=args.output_metrics)
