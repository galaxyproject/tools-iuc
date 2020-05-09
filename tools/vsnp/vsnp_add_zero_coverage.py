#!/usr/bin/env python

import argparse
import multiprocessing
import os
import pandas
import pysam
import queue
import re
import shutil
from Bio import SeqIO

INPUT_BAM_DIR = 'input_bam_dir'
INPUT_VCF_DIR = 'input_vcf_dir'
OUTPUT_VCF_DIR = 'output_vcf_dir'
OUTPUT_METRICS_DIR = 'output_metrics_dir'


def get_base_file_name(file_path):
    base_file_name = os.path.basename(file_path)
    if base_file_name.find(".") > 0:
        # Eliminate the extension.
        return os.path.splitext(base_file_name)[0]
    elif base_file_name.endswith("_vcf"):
        # The "." character has likely
        # changed to an "_" character.
        return base_file_name.rstrip("_vcf")
    return base_file_name


def get_coverage_and_snp_count(task_queue, reference, output_metrics, output_vcf, timeout):
    while True:
        try:
            tup = task_queue.get(block=True, timeout=timeout)
        except queue.Empty:
            break
        bam_file, vcf_file = tup
        # Create a coverage dictionary.
        coverage_dict = {}
        coverage_list = pysam.depth(bam_file, split_lines=True)
        for line in coverage_list:
            chrom, position, depth = line.split('\t')
            coverage_dict["%s-%s" % (chrom, position)] = depth
        # Convert it to a data frame.
        coverage_df = pandas.DataFrame.from_dict(coverage_dict, orient='index', columns=["depth"])
        # Create a zero coverage dictionary.
        zero_dict = {}
        for record in SeqIO.parse(reference, "fasta"):
            chrom = record.id
            total_len = len(record.seq)
            for pos in list(range(1, total_len + 1)):
                zero_dict["%s-%s" % (str(chrom), str(pos))] = 0
        # Convert it to a data frame with depth_x
        # and depth_y columns - index is NaN.
        zero_df = pandas.DataFrame.from_dict(zero_dict, orient='index', columns=["depth"])
        coverage_df = zero_df.merge(coverage_df, left_index=True, right_index=True, how='outer')
        # depth_x "0" column no longer needed.
        coverage_df = coverage_df.drop(columns=['depth_x'])
        coverage_df = coverage_df.rename(columns={'depth_y': 'depth'})
        # Covert the NaN to 0 coverage and get some metrics.
        coverage_df = coverage_df.fillna(0)
        coverage_df['depth'] = coverage_df['depth'].apply(int)
        total_length = len(coverage_df)
        average_coverage = coverage_df['depth'].mean()
        zero_df = coverage_df[coverage_df['depth'] == 0]
        total_zero_coverage = len(zero_df)
        total_coverage = total_length - total_zero_coverage
        genome_coverage = "{:.2%}".format(total_coverage / total_length)
        # Process the associated VCF input.
        column_names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample"]
        vcf_df = pandas.read_csv(vcf_file, sep='\t', header=None, names=column_names, comment='#')
        good_snp_count = len(vcf_df[(vcf_df['ALT'].str.len() == 1) & (vcf_df['REF'].str.len() == 1) & (vcf_df['QUAL'] > 150)])
        base_file_name = get_base_file_name(vcf_file)
        if total_zero_coverage > 0:
            header_file = "%s_header.csv" % base_file_name
            with open(header_file, 'w') as outfile:
                with open(vcf_file) as infile:
                    for line in infile:
                        if re.search('^#', line):
                            outfile.write("%s" % line)
            vcf_df_snp = vcf_df[vcf_df['REF'].str.len() == 1]
            vcf_df_snp = vcf_df_snp[vcf_df_snp['ALT'].str.len() == 1]
            vcf_df_snp['ABS_VALUE'] = vcf_df_snp['CHROM'].map(str) + "-" + vcf_df_snp['POS'].map(str)
            vcf_df_snp = vcf_df_snp.set_index('ABS_VALUE')
            cat_df = pandas.concat([vcf_df_snp, zero_df], axis=1, sort=False)
            cat_df = cat_df.drop(columns=['CHROM', 'POS', 'depth'])
            cat_df[['ID', 'ALT', 'QUAL', 'FILTER', 'INFO']] = cat_df[['ID', 'ALT', 'QUAL', 'FILTER', 'INFO']].fillna('.')
            cat_df['REF'] = cat_df['REF'].fillna('N')
            cat_df['FORMAT'] = cat_df['FORMAT'].fillna('GT')
            cat_df['Sample'] = cat_df['Sample'].fillna('./.')
            cat_df['temp'] = cat_df.index.str.rsplit('-', n=1)
            cat_df[['CHROM', 'POS']] = pandas.DataFrame(cat_df.temp.values.tolist(), index=cat_df.index)
            cat_df = cat_df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'Sample']]
            cat_df['POS'] = cat_df['POS'].astype(int)
            cat_df = cat_df.sort_values(['CHROM', 'POS'])
            body_file = "%s_body.csv" % base_file_name
            cat_df.to_csv(body_file, sep='\t', header=False, index=False)
            if output_vcf is None:
                output_vcf_file = os.path.join(OUTPUT_VCF_DIR, "%s.vcf" % base_file_name)
            else:
                output_vcf_file = output_vcf
            with open(output_vcf_file, "w") as outfile:
                for cf in [header_file, body_file]:
                    with open(cf, "r") as infile:
                        for line in infile:
                            outfile.write("%s" % line)
        else:
            if output_vcf is None:
                output_vcf_file = os.path.join(OUTPUT_VCF_DIR, "%s.vcf" % base_file_name)
            else:
                output_vcf_file = output_vcf
            shutil.copyfile(vcf_file, output_vcf_file)
        bam_metrics = [base_file_name, "", "%4f" % average_coverage, genome_coverage]
        vcf_metrics = [base_file_name, str(good_snp_count), "", ""]
        if output_metrics is None:
            output_metrics_file = os.path.join(OUTPUT_METRICS_DIR, "%s.tabular" % base_file_name)
        else:
            output_metrics_file = output_metrics
        metrics_columns = ["File", "Number of Good SNPs", "Average Coverage", "Genome Coverage"]
        with open(output_metrics_file, "w") as fh:
            fh.write("# %s\n" % "\t".join(metrics_columns))
            fh.write("%s\n" % "\t".join(bam_metrics))
            fh.write("%s\n" % "\t".join(vcf_metrics))
        task_queue.task_done()


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

    parser.add_argument('--output_metrics', action='store', dest='output_metrics', required=False, default=None, help='Output metrics text file')
    parser.add_argument('--output_vcf', action='store', dest='output_vcf', required=False, default=None, help='Output VCF file')
    parser.add_argument('--reference', action='store', dest='reference', help='Reference dataset')
    parser.add_argument('--processes', action='store', dest='processes', type=int, help='User-selected number of processes to use for job splitting')

    args = parser.parse_args()

    # The assumption here is that the list of files
    # in both INPUT_BAM_DIR and INPUT_VCF_DIR are
    # equal in number and named such that they are
    # properly matched if the directories contain
    # more than 1 file (i.e., hopefully the bam file
    # names and vcf file names will be something like
    # Mbovis-01D6_* so they can be # sorted and properly
    # associated with each other).
    bam_files = []
    for file_name in sorted(os.listdir(INPUT_BAM_DIR)):
        file_path = os.path.abspath(os.path.join(INPUT_BAM_DIR, file_name))
        bam_files.append(file_path)
    vcf_files = []
    for file_name in sorted(os.listdir(INPUT_VCF_DIR)):
        file_path = os.path.abspath(os.path.join(INPUT_VCF_DIR, file_name))
        vcf_files.append(file_path)

    multiprocessing.set_start_method('spawn')
    queue1 = multiprocessing.JoinableQueue()
    num_files = len(bam_files)
    cpus = set_num_cpus(num_files, args.processes)
    # Set a timeout for get()s in the queue.
    timeout = 0.05

    # Add each associated bam and vcf file pair to the queue.
    for i, bam_file in enumerate(bam_files):
        vcf_file = vcf_files[i]
        queue1.put((bam_file, vcf_file))

    # Complete the get_coverage_and_snp_count task.
    processes = [multiprocessing.Process(target=get_coverage_and_snp_count, args=(queue1, args.reference, args.output_metrics, args.output_vcf, timeout, )) for _ in range(cpus)]
    for p in processes:
        p.start()
    for p in processes:
        p.join()
    queue1.join()

    if queue1.empty():
        queue1.close()
        queue1.join_thread()
