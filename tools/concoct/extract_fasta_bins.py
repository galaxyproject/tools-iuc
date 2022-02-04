#!/usr/bin/env python

import argparse
import os
import sys

import pandas as pd
from Bio import SeqIO
from collections import defaultdict

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument("--input_fasta", action="store", dest="input_fasta", help="Input Fasta file")
parser.add_argument("--input_cluster", action="store", dest="input_cluster", help="Concoct output cluster file")
parser.add_argument("--output_path", help="Output directory")

args = parser.parse_args()

all_seqs = {}
for i, seq in enumerate(SeqIO.parse(args.input_fasta, "fasta")):
    all_seqs[seq.id] = seq
df = pd.read_csv(args.input_cluster)
try:
    assert df.columns[0] == 'contig_id'
    assert df.columns[1] == 'cluster_id'
except AssertionError:
    sys.stderr.write("ERROR! Header line was not 'contig_id, cluster_id', please adjust your input file. Exiting!\n")
    sys.exit(-1)

cluster_to_contigs = defaultdict(list)
for i, row in df.iterrows():
    cluster_to_contigs[row['cluster_id']].append(row['contig_id'])

for cluster_id, contig_ids in cluster_to_contigs.items():
    output_file = os.path.join(args.output_path, "{0}.fa".format(cluster_id))
    seqs = [all_seqs[contig_id] for contig_id in contig_ids]
    with open(output_file, 'w') as ofh:
        SeqIO.write(seqs, ofh, 'fasta')
