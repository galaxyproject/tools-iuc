#!/usr/bin/env python
"""
Input: DNA FASTA file + Glimmer ORF file
Output: ORF sequences as FASTA file
Author: Bjoern Gruening
"""
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def glimmer2seq(glimmer_prediction=sys.argv[1], genome_sequence=sys.argv[2], outfile=sys.argv[3]):
    if len(sys.argv) >= 4:
        glimmerfile = open(glimmer_prediction, "r")
        sequence = open(genome_sequence)
    else:
        print("Missing input values.")
        sys.exit()

    fastafile = SeqIO.parse(sequence, "fasta")

    sequences = dict()
    seq_records = list()
    for entry in fastafile:
        sequences[entry.description] = entry

    for line in glimmerfile:
        if line.startswith('>'):
            entry = sequences[line[1:].strip()]
        else:
            orf_start = int(line[8:17])
            orf_end = int(line[18:26])

            orf_name = line[0:8]
            if orf_start <= orf_end:
                seq_records.append(SeqRecord(entry.seq[orf_start - 1:orf_end], id=orf_name, description=entry.description))
            else:
                seq_records.append(SeqRecord(entry.seq[orf_end - 1:orf_start].reverse_complement(), id=orf_name, description=entry.description))

    SeqIO.write(seq_records, outfile, "fasta")
    glimmerfile.close()
    sequence.close()


if __name__ == "__main__":
    glimmer2seq()
