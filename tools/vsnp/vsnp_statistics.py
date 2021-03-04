#!/usr/bin/env python

import argparse
import csv
import gzip
import os
from functools import partial

import numpy
import pandas
from Bio import SeqIO


def nice_size(size):
    # Returns a readably formatted string with the size
    words = ['bytes', 'KB', 'MB', 'GB', 'TB', 'PB', 'EB']
    prefix = ''
    try:
        size = float(size)
        if size < 0:
            size = abs(size)
            prefix = '-'
    except Exception:
        return '??? bytes'
    for ind, word in enumerate(words):
        step = 1024 ** (ind + 1)
        if step > size:
            size = size / float(1024 ** ind)
            if word == 'bytes':  # No decimals for bytes
                return "%s%d bytes" % (prefix, size)
            return "%s%.1f %s" % (prefix, size, word)
    return '??? bytes'


def output_statistics(fastq_files, idxstats_files, metrics_files, output_file, gzipped, dbkey):
    # Produce an Excel spreadsheet that
    # contains a row for each sample.
    columns = ['Reference', 'File Size', 'Mean Read Length', 'Mean Read Quality', 'Reads Passing Q30',
               'Total Reads', 'All Mapped Reads', 'Unmapped Reads', 'Unmapped Reads Percentage of Total',
               'Reference with Coverage', 'Average Depth of Coverage', 'Good SNP Count']
    data_frames = []
    for i, fastq_file in enumerate(fastq_files):
        idxstats_file = idxstats_files[i]
        metrics_file = metrics_files[i]
        file_name_base = os.path.basename(fastq_file)
        # Read fastq_file into a data frame.
        _open = partial(gzip.open, mode='rt') if gzipped else open
        with _open(fastq_file) as fh:
            identifiers = []
            seqs = []
            letter_annotations = []
            for seq_record in SeqIO.parse(fh, "fastq"):
                identifiers.append(seq_record.id)
                seqs.append(seq_record.seq)
                letter_annotations.append(seq_record.letter_annotations["phred_quality"])
        # Convert lists to Pandas series.
        s1 = pandas.Series(identifiers, name='id')
        s2 = pandas.Series(seqs, name='seq')
        # Gather Series into a data frame.
        fastq_df = pandas.DataFrame(dict(id=s1, seq=s2)).set_index(['id'])
        total_reads = int(len(fastq_df.index) / 4)
        current_sample_df = pandas.DataFrame(index=[file_name_base], columns=columns)
        # Reference
        current_sample_df.at[file_name_base, 'Reference'] = dbkey
        # File Size
        current_sample_df.at[file_name_base, 'File Size'] = nice_size(os.path.getsize(fastq_file))
        # Mean Read Length
        sampling_size = 10000
        if sampling_size > total_reads:
            sampling_size = total_reads
        fastq_df = fastq_df.iloc[3::4].sample(sampling_size)
        dict_mean = {}
        list_length = []
        i = 0
        for id, seq, in fastq_df.iterrows():
            dict_mean[id] = numpy.mean(letter_annotations[i])
            list_length.append(len(seq.array[0]))
            i += 1
        current_sample_df.at[file_name_base, 'Mean Read Length'] = '%.1f' % numpy.mean(list_length)
        # Mean Read Quality
        df_mean = pandas.DataFrame.from_dict(dict_mean, orient='index', columns=['ave'])
        current_sample_df.at[file_name_base, 'Mean Read Quality'] = '%.1f' % df_mean['ave'].mean()
        # Reads Passing Q30
        reads_gt_q30 = len(df_mean[df_mean['ave'] >= 30])
        reads_passing_q30 = '{:10.2f}'.format(reads_gt_q30 / sampling_size)
        current_sample_df.at[file_name_base, 'Reads Passing Q30'] = reads_passing_q30
        # Total Reads
        current_sample_df.at[file_name_base, 'Total Reads'] = total_reads
        # All Mapped Reads
        all_mapped_reads, unmapped_reads = process_idxstats_file(idxstats_file)
        current_sample_df.at[file_name_base, 'All Mapped Reads'] = all_mapped_reads
        # Unmapped Reads
        current_sample_df.at[file_name_base, 'Unmapped Reads'] = unmapped_reads
        # Unmapped Reads Percentage of Total
        if unmapped_reads > 0:
            unmapped_reads_percentage = '{:10.2f}'.format(unmapped_reads / total_reads)
        else:
            unmapped_reads_percentage = 0
        current_sample_df.at[file_name_base, 'Unmapped Reads Percentage of Total'] = unmapped_reads_percentage
        # Reference with Coverage
        ref_with_coverage, avg_depth_of_coverage, good_snp_count = process_metrics_file(metrics_file)
        current_sample_df.at[file_name_base, 'Reference with Coverage'] = ref_with_coverage
        # Average Depth of Coverage
        current_sample_df.at[file_name_base, 'Average Depth of Coverage'] = avg_depth_of_coverage
        # Good SNP Count
        current_sample_df.at[file_name_base, 'Good SNP Count'] = good_snp_count
        data_frames.append(current_sample_df)
    output_df = pandas.concat(data_frames)
    output_df.to_csv(output_file, sep='\t', quoting=csv.QUOTE_NONE, escapechar='\\')


def process_idxstats_file(idxstats_file):
    all_mapped_reads = 0
    unmapped_reads = 0
    with open(idxstats_file, "r") as fh:
        for i, line in enumerate(fh):
            line = line.rstrip('\r\n')
            items = line.split("\t")
            if i == 0:
                # NC_002945.4 4349904 213570 4047
                all_mapped_reads = int(items[2])
            elif i == 1:
                # * 0 0 82774
                unmapped_reads = int(items[3])
    return all_mapped_reads, unmapped_reads


def process_metrics_file(metrics_file):
    ref_with_coverage = '0%'
    avg_depth_of_coverage = 0
    good_snp_count = 0
    with open(metrics_file, "r") as ifh:
        for i, line in enumerate(ifh):
            if i == 0:
                # Skip comments.
                continue
            line = line.rstrip('\r\n')
            items = line.split("\t")
            if i == 1:
                # MarkDuplicates 10.338671 98.74%
                ref_with_coverage = items[3]
                avg_depth_of_coverage = items[2]
            elif i == 2:
                # VCFfilter 611
                good_snp_count = items[1]
    return ref_with_coverage, avg_depth_of_coverage, good_snp_count


parser = argparse.ArgumentParser()

parser.add_argument('--dbkey', action='store', dest='dbkey', help='Reference dbkey')
parser.add_argument('--gzipped', action='store_true', dest='gzipped', required=False, default=False, help='Input files are gzipped')
parser.add_argument('--input_idxstats_dir', action='store', dest='input_idxstats_dir', required=False, default=None, help='Samtools idxstats input directory')
parser.add_argument('--input_metrics_dir', action='store', dest='input_metrics_dir', required=False, default=None, help='vSNP add zero coverage metrics input directory')
parser.add_argument('--input_reads_dir', action='store', dest='input_reads_dir', required=False, default=None, help='Samples input directory')
parser.add_argument('--list_paired', action='store_true', dest='list_paired', required=False, default=False, help='Input samples is a list of paired reads')
parser.add_argument('--output', action='store', dest='output', help='Output Excel statistics file')
parser.add_argument('--read1', action='store', dest='read1', help='Required: single read')
parser.add_argument('--read2', action='store', dest='read2', required=False, default=None, help='Optional: paired read')
parser.add_argument('--samtools_idxstats', action='store', dest='samtools_idxstats', help='Output of samtools_idxstats')
parser.add_argument('--vsnp_azc', action='store', dest='vsnp_azc', help='Output of vsnp_add_zero_coverage')

args = parser.parse_args()

fastq_files = []
idxstats_files = []
metrics_files = []
# Accumulate inputs.
if args.read1 is not None:
    # The inputs are not dataset collections, so
    # read1, read2 (possibly) and vsnp_azc will also
    # not be None.
    fastq_files.append(args.read1)
    idxstats_files.append(args.samtools_idxstats)
    metrics_files.append(args.vsnp_azc)
    if args.read2 is not None:
        fastq_files.append(args.read2)
        idxstats_files.append(args.samtools_idxstats)
        metrics_files.append(args.vsnp_azc)
else:
    for file_name in sorted(os.listdir(args.input_reads_dir)):
        fastq_files.append(os.path.join(args.input_reads_dir, file_name))
    for file_name in sorted(os.listdir(args.input_idxstats_dir)):
        idxstats_files.append(os.path.join(args.input_idxstats_dir, file_name))
        if args.list_paired:
            # Add the idxstats file for reverse.
            idxstats_files.append(os.path.join(args.input_idxstats_dir, file_name))
    for file_name in sorted(os.listdir(args.input_metrics_dir)):
        metrics_files.append(os.path.join(args.input_metrics_dir, file_name))
        if args.list_paired:
            # Add the metrics file for reverse.
            metrics_files.append(os.path.join(args.input_metrics_dir, file_name))
output_statistics(fastq_files, idxstats_files, metrics_files, args.output, args.gzipped, args.dbkey)
