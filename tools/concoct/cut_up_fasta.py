#!/usr/bin/env python

import argparse
import gzip
from functools import partial

from Bio import SeqIO


def cut_up_fasta(input_fasta, chunk_size, overlap, merge_last, output_fasta, output_bed, gzipped):
    if gzipped:
        _open = partial(gzip.open, mode='rt')
    else:
        _open = open

    fasta_fh = open(output_fasta, 'w')

    if output_bed is not None:
        bed_fh = open(output_bed, 'w')

    with _open(input_fasta) as input_fh:
        for record in SeqIO.parse(input_fh, "fasta"):
            if (not merge_last and len(record.seq) > chunk_size) or (merge_last and len(record.seq) >= 2 * chunk_size):
                for index, split_seq in enumerate(chunks(record.seq, chunk_size, overlap, merge_last)):
                    fasta_fh.write(f">{record.id}.{index}\n{split_seq}\n")
                    if output_bed is not None:
                        bed_fh.write(f"{record.id}\t{chunk_size * index}\t{chunk_size * index + len(split_seq)}\t{record.id}.{index}\n")
            else:
                fasta_fh.write(f">{record.id}\n{record.seq}\n")
                if output_bed is not None:
                    bed_fh.write(f"{record.id}\t0\t{len(record.seq)}\t{record.id}\n")
    if output_bed is not None:
        bed_fh.close()


def chunks(seq, chunk_size, overlap_size, merge_last):
    # Yield successive chunk_size-sized chunks from seq
    # with given overlap overlap_size between the chunks.
    assert chunk_size > overlap_size
    if merge_last:
        for i in range(0, len(seq) - chunk_size + 1, chunk_size - overlap_size):
            yield seq[i:i + chunk_size] if i + chunk_size + chunk_size - overlap_size <= len(seq) else seq[i:]
    else:
        for i in range(0, len(seq), chunk_size - overlap_size):
            yield seq[i:i + chunk_size]


parser = argparse.ArgumentParser()
parser.add_argument("--input_fasta", action="store", dest="input_fasta", help="Fasta files with contigs")
parser.add_argument("--gzipped", action="store_true", dest="gzipped", help="Input file is gzipped")
parser.add_argument("--chunk_size", action="store", dest="chunk_size", type=int, help="Chunk size\n")
parser.add_argument("--overlap_size", action="store", dest="overlap_size", type=int, help="Overlap size\n")
parser.add_argument("--merge_last", default=False, action="store_true", dest="merge_last", help="Concatenate final part to last contig\n")
parser.add_argument("--output_bed", action="store", dest="output_bed", default=None, help="BED file to be created with exact regions of the original contigs corresponding to the newly created contigs")
parser.add_argument("--output_fasta", action="store", dest="output_fasta", help="Output fasta file with cut contigs")

args = parser.parse_args()
cut_up_fasta(args.input_fasta, args.chunk_size, args.overlap_size, args.merge_last, args.output_fasta, args.output_bed, args.gzipped)
