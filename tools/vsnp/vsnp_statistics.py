#!/usr/bin/env python

import argparse
import gzip
import os
from functools import partial

import numpy
import pandas
from Bio import SeqIO


class Statistics:

    def __init__(self, reference, fastq_file, file_size, total_reads, mean_read_length, mean_read_quality, reads_passing_q30):
        self.reference = reference
        self.fastq_file = fastq_file
        self.file_size = file_size
        self.total_reads = total_reads
        self.mean_read_length = mean_read_length
        self.mean_read_quality = mean_read_quality
        self.reads_passing_q30 = reads_passing_q30


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


def get_statistics(dbkey, fastq_file, gzipped):
    sampling_size = 10000
    # Read fastq_file into a data fram to
    # get the phred quality scores.
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
    # Starting at row 3, keep every 4 row
    # random sample specified number of rows.
    file_size = nice_size(os.path.getsize(fastq_file))
    total_reads = len(seqs)
    # Mean Read Length
    if sampling_size > total_reads:
        sampling_size = total_reads
    try:
        fastq_df = fastq_df.iloc[3::4].sample(sampling_size)
    except ValueError:
        fastq_df = fastq_df.iloc[3::4].sample(sampling_size, replace=True)
    dict_mean = {}
    list_length = []
    i = 0
    for id, seq, in fastq_df.iterrows():
        dict_mean[id] = numpy.mean(letter_annotations[i])
        list_length.append(len(seq.array[0]))
        i += 1
    mean_read_length = '%.1f' % numpy.mean(list_length)
    # Mean Read Quality
    df_mean = pandas.DataFrame.from_dict(dict_mean, orient='index', columns=['ave'])
    mean_read_quality = '%.1f' % df_mean['ave'].mean()
    # Reads Passing Q30
    reads_gt_q30 = len(df_mean[df_mean['ave'] >= 30])
    reads_passing_q30 = '{:10.2f}'.format(reads_gt_q30 / sampling_size)
    stats = Statistics(dbkey, os.path.basename(fastq_file), file_size, total_reads, mean_read_length,
                       mean_read_quality, reads_passing_q30)
    return stats


def accrue_statistics(dbkey, read1, read2, gzipped):
    read1_stats = get_statistics(dbkey, read1, gzipped)
    if read2 is None:
        read2_stats = None
    else:
        read2_stats = get_statistics(dbkey, read2, gzipped)
    return read1_stats, read2_stats


def output_statistics(read1_stats, read2_stats, idxstats_file, metrics_file, output_file):
    paired_reads = read2_stats is not None
    if paired_reads:
        columns = ['Read1 FASTQ', 'File Size', 'Reads', 'Mean Read Length', 'Mean Read Quality',
                   'Reads Passing Q30', 'Read2 FASTQ', 'File Size', 'Reads', 'Mean Read Length', 'Mean Read Quality',
                   'Reads Passing Q30', 'Total Reads', 'All Mapped Reads', 'Unmapped Reads',
                   'Unmapped Reads Percentage of Total', 'Reference with Coverage', 'Average Depth of Coverage',
                   'Good SNP Count', 'Reference']
    else:
        columns = ['FASTQ', 'File Size', 'Mean Read Length', 'Mean Read Quality', 'Reads Passing Q30',
                   'Total Reads', 'All Mapped Reads', 'Unmapped Reads', 'Unmapped Reads Percentage of Total',
                   'Reference with Coverage', 'Average Depth of Coverage', 'Good SNP Count', 'Reference']
    with open(output_file, "w") as outfh:
        # Make sure the header starts with a # so
        # MultiQC can properly handle the output.
        outfh.write("%s\n" % "\t".join(columns))
        line_items = []
        # Get the current stats and associated files.
        # Get and output the statistics.
        line_items.append(read1_stats.fastq_file)
        line_items.append(read1_stats.file_size)
        if paired_reads:
            line_items.append(read1_stats.total_reads)
        line_items.append(read1_stats.mean_read_length)
        line_items.append(read1_stats.mean_read_quality)
        line_items.append(read1_stats.reads_passing_q30)
        if paired_reads:
            line_items.append(read2_stats.fastq_file)
            line_items.append(read2_stats.file_size)
            line_items.append(read2_stats.total_reads)
            line_items.append(read2_stats.mean_read_length)
            line_items.append(read2_stats.mean_read_quality)
            line_items.append(read2_stats.reads_passing_q30)
        # Total Reads
        if paired_reads:
            total_reads = read1_stats.total_reads + read2_stats.total_reads
        else:
            total_reads = read1_stats.total_reads
        line_items.append(total_reads)
        # All Mapped Reads
        all_mapped_reads, unmapped_reads = process_idxstats_file(idxstats_file)
        line_items.append(all_mapped_reads)
        line_items.append(unmapped_reads)
        # Unmapped Reads Percentage of Total
        if unmapped_reads > 0:
            unmapped_reads_percentage = '{:10.2f}'.format(unmapped_reads / total_reads)
        else:
            unmapped_reads_percentage = 0
        line_items.append(unmapped_reads_percentage)
        # Reference with Coverage
        ref_with_coverage, avg_depth_of_coverage, good_snp_count = process_metrics_file(metrics_file)
        line_items.append(ref_with_coverage)
        line_items.append(avg_depth_of_coverage)
        line_items.append(good_snp_count)
        line_items.append(read1_stats.reference)
        outfh.write('%s\n' % '\t'.join(str(x) for x in line_items))


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
parser.add_argument('--output', action='store', dest='output', help='Output Excel statistics file')
parser.add_argument('--read1', action='store', dest='read1', help='Required: single read')
parser.add_argument('--read2', action='store', dest='read2', required=False, default=None, help='Optional: paired read')
parser.add_argument('--samtools_idxstats', action='store', dest='samtools_idxstats', help='Output of samtools_idxstats')
parser.add_argument('--vsnp_azc_metrics', action='store', dest='vsnp_azc_metrics', help='Output of vsnp_add_zero_coverage')

args = parser.parse_args()

stats_list = []
idxstats_files = []
metrics_files = []
# Accumulate inputs.
read1_stats, read2_stats = accrue_statistics(args.dbkey, args.read1, args.read2, args.gzipped)
output_statistics(read1_stats, read2_stats, args.samtools_idxstats, args.vsnp_azc_metrics, args.output)
