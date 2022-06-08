#!/usr/bin/env python

import argparse
import os


class Statistics:

    def __init__(self, file_name, file_size, seq_type, num_seqs, sum_len, min_len, avg_len,
                 max_len, q1, q2, q3, sum_gap, n50, pass_q20, pass_q30, read_quality_average):
        self.file_name = file_name
        self.file_size = file_size
        self.seq_type = seq_type
        self.num_seqs = num_seqs
        self.sum_len = sum_len
        self.min_len = min_len
        self.avg_len = avg_len
        self.max_len = max_len
        self.q1 = q1
        self.q2 = q2
        self.q3 = q3
        self.sum_gap = sum_gap
        self.n50 = n50
        self.pass_q20 = pass_q20
        self.pass_q30 = pass_q30
        self.read_quality_average = read_quality_average


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


def output_statistics(read1_stats, read2_stats, output_file):
    paired_reads = read2_stats is not None
    if paired_reads:
        columns = ['R1 FASTQ', 'R1 File Size', 'R1 Read Count', 'R1 Length Sum', 'R1 Min Length',
                   'R1 Ave Length', 'R1 Max Length', 'R1 Q1', 'R1 Q2', 'R1 Q3', 'R1 Sum Gap',
                   'R1 N50', 'R1 Passing Q20', 'R1 Passing Q30', 'R1 Read Quality Ave', 'R2 FASTQ',
                   'R2 File Size', 'R2 Read Count', 'R2 Length Sum', 'R2 Min Length', 'R2 Ave Length',
                   'R2 Max Length', 'R2 Q1', 'R2 Q2', 'R2 Q3', 'R2 Sum Gap', 'R2 N50', 'R2 Passing Q20',
                   'R2 Passing Q30', 'R2 Read Quality Ave']
    else:
        columns = ['FASTQ', 'File Size', 'Read Count', 'Length Sum', 'Min Length', 'Ave Length',
                   'Max Length', 'Q1', 'Q2', 'Q3', 'Sum Gap', 'N50', 'Passing Q20', 'Passing Q30',
                   'Read Quality Ave']
    with open(output_file, "w") as outfh:
        # Make sure the header starts with a # so
        # MultiQC can properly handle the output.
        outfh.write("%s\n" % "\t".join(columns))
        line_items = []
        # Get the current stats and associated files.
        # Get and output the statistics.
        line_items.append(read1_stats.file_name)
        line_items.append(read1_stats.file_size)
        line_items.append(read1_stats.num_seqs)
        line_items.append(read1_stats.sum_len)
        line_items.append(read1_stats.min_len)
        line_items.append(read1_stats.avg_len)
        line_items.append(read1_stats.max_len)
        line_items.append(read1_stats.q1)
        line_items.append(read1_stats.q2)
        line_items.append(read1_stats.q3)
        line_items.append(read1_stats.sum_gap)
        line_items.append(read1_stats.n50)
        line_items.append(read1_stats.pass_q20)
        line_items.append(read1_stats.pass_q30)
        line_items.append(read1_stats.read_quality_average)
        if paired_reads:
            line_items.append(read2_stats.file_name)
            line_items.append(read2_stats.file_size)
            line_items.append(read2_stats.num_seqs)
            line_items.append(read2_stats.sum_len)
            line_items.append(read2_stats.min_len)
            line_items.append(read2_stats.avg_len)
            line_items.append(read2_stats.max_len)
            line_items.append(read2_stats.q1)
            line_items.append(read2_stats.q2)
            line_items.append(read2_stats.q3)
            line_items.append(read2_stats.sum_gap)
            line_items.append(read2_stats.n50)
            line_items.append(read2_stats.pass_q20)
            line_items.append(read2_stats.pass_q30)
            line_items.append(read2_stats.read_quality_average)
        outfh.write('%s\n' % '\t'.join(str(x) for x in line_items))


def get_statistics(fastq_file, seqkit_stats_file, seqkit_fx2tab_file):
    file_size = nice_size(os.path.getsize(fastq_file))
    # SeqKit statistics.
    with open(seqkit_stats_file, "r") as fh:
        # This is a 2-line file
        for i, line in enumerate(fh):
            if i == 0:
                # Skip header
                continue
            line = line.rstrip('\r\n')
            items = line.split("\t")
            file_name = fastq_file
            seq_type = items[2]
            num_seqs = items[3]
            sum_len = items[4]
            min_len = items[5]
            avg_len = items[6]
            max_len = items[7]
            q1 = items[8]
            q2 = items[9]
            q3 = items[10]
            sum_gap = items[11]
            n50 = items[12]
            try:
                pass_q20 = items[13]
            except IndexError:
                pass_q20 = 0
            try:
                pass_q30 = items[14]
            except IndexError:
                pass_q30 = 0
    # Average read quality is not normalized on length.
    avg_sum = 0
    with open(seqkit_fx2tab_file, "r") as fh:
        for i, line in enumerate(fh):
            if i == 0:
                # Skip header
                continue
            line = line.rstrip('\r\n')
            items = line.split("\t")
            avg_sum += float(items[3])
    read_quality_average = "{:.2f}".format(avg_sum / float(i - 1))
    return Statistics(file_name, file_size, seq_type, num_seqs, sum_len, min_len, avg_len,
                      max_len, q1, q2, q3, sum_gap, n50, pass_q20, pass_q30, read_quality_average)


parser = argparse.ArgumentParser()

parser.add_argument('--output', action='store', dest='output', help='Output Excel statistics file')
parser.add_argument('--read1', action='store', dest='read1', help='Required: single read')
parser.add_argument('--read2', action='store', dest='read2', required=False, default=None, help='Optional: paired read')
parser.add_argument('--read1_seqkit_stats', action='store', dest='read1_seqkit_stats', help='Output of SeqKit statistics for forward read')
parser.add_argument('--read2_seqkit_stats', action='store', dest='read2_seqkit_stats', required=False, default=None, help='Output of SeqKit statistics for reverse read')
parser.add_argument('--read1_seqkit_fx2tab', action='store', dest='read1_seqkit_fx2tab', help='Output of SeqKit fx2tab for forward read')
parser.add_argument('--read2_seqkit_fx2tab', action='store', dest='read2_seqkit_fx2tab', required=False, default=None, help='Output of SeqKit fx2tab for reverse read')

args = parser.parse_args()

read1_stats = get_statistics(args.read1, args.read1_seqkit_stats, args.read1_seqkit_fx2tab)
if args.read2 is None:
    read2_stats = None
else:
    read2_stats = get_statistics(args.read2, args.read2_seqkit_stats, args.read2_seqkit_fx2tab)

output_statistics(read1_stats, read2_stats, args.output)
