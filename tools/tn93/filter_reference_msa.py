import csv
import argparse
from Bio import SeqIO

# python filter_reference_msa.py --reference '$reference_msa' --distances pairwise.csv --output '$filtered_reference'


arguments = argparse.ArgumentParser(description='Combine alignments into a single file, adding a reference sequence as well')

arguments.add_argument('-f', '--reference', help='Reference sequence', required=True, type=str)
arguments.add_argument('-d', '--distances', help='Calculated pairwise distances', required=True, type=str)
arguments.add_argument('-r', '--reads', help='Output file for filtered reads', required=True, type=str)
arguments.add_argument('-q', '--clusters', help='Compressed clusters', required=True, type=str)
settings = arguments.parse_args()

reference_name = 'REFERENCE'
reference_seq = ''

with open(settings.reference) as seq_fh:
    for seq_record in SeqIO.parse(seq_fh, 'fasta'):
        reference_name = seq_record.name
        reference_seq = seq_record.seq
        print('Found reference name: %s' % reference_name)
        break

with open(settings.distances) as fh:
    reader = csv.reader(fh, delimiter=',')
    next(reader)
    seqs_to_filter = set()
    for line in reader:
        if line[1] not in seqs_to_filter:
            seqs_to_filter.add(line[1])
    if reference_name in seqs_to_filter:
        seqs_to_filter.remove(reference_name)

with open(settings.reads, "a+") as fh:
    seqs_filtered = list()
    for seq_record in SeqIO.parse(settings.clusters, "fasta"):
        if seq_record.name not in seqs_to_filter:
            if seq_record.name == reference_name:
                if seq_record.name not in seqs_filtered:
                    seqs_filtered.append(seq_record.name)
                else:
                    continue
                print('Replacing %s with REFERENCE because it matches %s' % (seq_record.name, reference_name))
                print("\n>%s\n%s" % ("REFERENCE", str(seq_record.seq)), file=fh)
            else:
                print('Not replacing %s with REFERENCE because it does not match %s' % (seq_record.name, reference_name))
                print("\n>%s\n%s" % (seq_record.name, str(seq_record.seq)), file=fh)
    if reference_name not in seqs_filtered:
        fh.write('\n>REFERENCE\n%s' % reference_seq)
