import argparse
import csv
import random

from Bio import SeqIO

arguments = argparse.ArgumentParser(description='Combine alignments into a single file, adding a reference sequence as well')

arguments.add_argument('-f', '--reference', help='Reference sequence', required=True, type=str)
arguments.add_argument('-d', '--distances', help='Calculated pairwise distances', required=True, type=str)
arguments.add_argument('-r', '--reads', help='Output file for filtered reads', required=True, type=str)
arguments.add_argument('-q', '--clusters', help='Compressed background clusters', required=True, type=str)
settings = arguments.parse_args()

reference_name = 'REFERENCE'
reference_seq = ''


def unique_id(new_id, existing_ids):
    while new_id in existing_ids:
        new_id += '_' + ''.join(random.choices('0123456789abcdef', k=10))
    return new_id


with open(settings.reference) as seq_fh:
    for seq_record in SeqIO.parse(seq_fh, 'fasta'):
        reference_name = seq_record.name.split(' ')[0]
        reference_seq = seq_record.seq
        break

with open(settings.distances) as fh:
    reader = csv.reader(fh, delimiter=',')
    next(reader)
    seqs_to_filter = set()
    for line in reader:
        if line[1] not in seqs_to_filter:
            seqs_to_filter.add(line[1])
        else:
            seqs_to_filter.add(unique_id(line[1], seqs_to_filter))
    if reference_name in seqs_to_filter:
        seqs_to_filter.remove(reference_name)

with open(settings.reads, "a+") as fh:
    seqs_filtered = list()
    for seq_record in SeqIO.parse(settings.clusters, "fasta"):
        if seq_record.name.split(' ')[0] == reference_name:
            continue
        if seq_record.name not in seqs_to_filter:
            unique_name = unique_id(seq_record.name, seqs_filtered)
            fh.write('\n>%s\n%s' % (unique_name, seq_record.seq))
            seqs_filtered.append(unique_name)
    if reference_name not in seqs_filtered:
        fh.write('\n>REFERENCE\n%s' % reference_seq)
