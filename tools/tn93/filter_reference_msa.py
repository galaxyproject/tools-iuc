import csv
import argparse
from Bio import SeqIO

# python filter_reference_msa.py --reference '$reference_msa' --distances pairwise.csv --output '$filtered_reference'


arguments = argparse.ArgumentParser(description='Combine alignments into a single file, adding a reference sequence as well')

arguments.add_argument('-r', '--reference', help='Reference alignment to filter', required=True, type=str )
arguments.add_argument('-d', '--distances', help='Calculated pairwise distances', required=True, type=str)
arguments.add_argument('-o', '--output', help='Filtered references', required=True, type=str)
arguments.add_argument('-s', '--sequence', help='Reference sequence', required=True, type=str)
settings = arguments.parse_args()

with open(options.sequence) as seq_fh:
    for seq_record in SeqIO.parse(seq_fh, 'fasta'):
        reference_name = seq_record.name
        break

with open (settings.distances) as fh:
    reader = csv.reader (fh, delimiter = ',')
    next (reader)
    seqs_to_filter = set ()
    for l in reader:
        seqs_to_filter.add (l[1])
    if reference_name in seqs_to_filter:
        seqs_to_filter.remove (reference_name)
with open (options.output, "a+") as fh:
    for seq_record in SeqIO.parse(options.reference, "fasta"):
        if not seq_record.name in seqs_to_filter:
            if seq_record.name == reference_name:
                print ("\n>%s\n%s" % ("REFERENCE", str(seq_record.seq)), file = fh)
            else:
                print ("\n>%s\n%s" % (seq_record.name, str(seq_record.seq)), file = fh)

