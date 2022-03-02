#!/usr/bin/env python

import argparse
import gzip
from functools import partial

from Bio import SeqIO


def cut_up_fasta(input_fasta, chunk_size, overlap, merge_last, bedoutfile, gzipped):
    if args.gzipped:
        _open = partial(gzip.open, mode='rt')
    else:
        _open = open

    if bedoutfile:
        bedoutfile_fh = open(bedoutfile, 'w')

    with _open(input_fasta) as fh:
        for record in SeqIO.parse(fh, "fasta"):
            if (not merge_last and len(record.seq) > chunk_size) or (merge_last and len(record.seq) >= 2 * chunk_size):
                i = 0
                for split_seq in chunks(record.seq, chunk_size, overlap, merge_last):
                    print(">%s.%i\n%s" % (record.id, i, split_seq))
                    if bedoutfile:
                        print("{0}\t{2}\t{3}\t{0}.{1}".format(record.id, i, chunk_size * i, chunk_size * i + len(split_seq)),
                              file=bedoutfile_fh)
                    i = i + 1
            else:
                print(">%s\n%s" % (record.id, record.seq))
                if bedoutfile:
                    print("{0}\t0\t{1}\t{0}".format(record.id, len(record.seq)),
                          file=bedoutfile_fh)

    if bedoutfile:
        bedoutfile_fh.close()


def chunks(seq, chunk_size, overlap_size, merge_last):
    # Yield successive chunk_size-sized chunks from l
    # with given overlap o between the chunks.
    assert chunk_size > overlap_size

    if not merge_last:
        for i in range(0, len(seq), chunk_size - overlap_size):
            yield seq[i:i + chunk_size]
    else:
        for i in range(0, len(seq) - chunk_size + 1, chunk_size - overlap_size):
            yield seq[i:i + chunk_size] if i + chunk_size + chunk_size - overlap_size <= len(seq) else seq[i:]


parser = argparse.ArgumentParser()
parser.add_argument("--input_fasta", action="store", dest="input_fasta", help="Fasta files with contigs")
parser.add_argument("--gzipped", action="store_true", dest="gzipped", help="Input file is gzipped")
parser.add_argument("--chunk_size", action="store", dest="chunk_size", type=int, help="Chunk size\n")
parser.add_argument("--overlap_size",action="store", dest="overlap_size", type=int, help="Overlap size\n")
parser.add_argument("--merge_last", default=False, action="store_true", dest="merge_last", help="Concatenate final part to last contig\n")
parser.add_argument("--bedfile", action="store", dest="bedfile", default=None, help="BEDfile to be created with exact regions of the original contigs corresponding to the newly created contigs")

args = parser.parse_args()
cut_up_fasta(args.input_fasta, args.chunk_size, args.overlap_size, args.merge_last, args.bedfile, args.gzipped)
