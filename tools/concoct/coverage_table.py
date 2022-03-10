#!/usr/bin/env python

import argparse
import gzip
from functools import partial

from Bio import SeqIO


def generate_coverage_table(input_fasta, input_tabular, gzipped, output):
    # Read input file into a dict and return everything
    # in the table format required by CONCOCT.
    gc_and_len_dict = get_gc_and_len_dict(input_fasta, gzipped)
    assert(len(gc_and_len_dict) > 0)
    bed_coverage_dict = get_bed_coverage_dict(input_tabular)

    with open(output, 'w') as fh:
        # Output the header.
        fh.write("contig\tlength")
        t = tuple(range(len(bed_coverage_dict)))
        fh.write("\tcov_mean_sample_%d\n" % len(t))
        # Output the content.
        for acc in gc_and_len_dict:
            # Fasta stats.
            fh.write("%s\t%s" % (acc, gc_and_len_dict[acc]['length']))
            # Mean
            try:
                # Coverage mean
                fh.write("\t%f" % (bed_coverage_dict[acc]["cov_mean"]))
            except KeyError:
                # No reads mapped to this contig
                fh.write("\t0")
            fh.write("\n")


def get_bed_coverage_dict(input_tabular):
    # Ddetermine mean coverage and percentage covered
    # for each contig, returning a dict with fasta id
    # as key and percentage covered and cov_mean as keys
    # for the inner dict.
    out_dict = {}
    debug_fh = open("/tmp/coverage.log", "w")

    with open(input_tabular, 'r') as fh:
        for line in fh:
            line = line.rstrip('\r\n')
            debug_fh.write("\nline:\n%s" % str(line))
            cols = line.split('\t')
            debug_fh.write("\ncols:\n%s" % str(cols))
            try:
                d = out_dict[cols[0]]
            except KeyError:
                d = {}
                out_dict[cols[0]] = d
            debug_fh.write("\nint(cols[1]): %d" % int(cols[1]))
            debug_fh.write("\nfloat(cols[4]): %d" % float(cols[4]))
            if int(cols[1]) == 0:
                d["percentage_covered"] = 100 - float(cols[4]) * 100.0
            else:
                d["cov_mean"] = d.get("cov_mean", 0) + int(cols[1]) * float(cols[4])
    debug_fh.close()
    return out_dict


def get_gc_and_len_dict(input_fasta, gzipped):
    # Creates a dictionary with the fasta id as key
    # and GC and length as keys for the inner dictionary.
    if gzipped:
        _open = partial(gzip.open, mode='rt')
    else:
        _open = open

    out_dict = {}
    with _open(input_fasta) as input_fh:
        for rec in SeqIO.parse(input_fh, "fasta"):
            out_dict[rec.id] = {}
            out_dict[rec.id]["length"] = len(rec.seq)
    return out_dict


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input_tabular', action='store', dest='input_tabular', help='bedtools genomeCoverageBed bed file')
parser.add_argument('--input_fasta', action='store', dest='input_fasta', help='Contigs fasta file')
parser.add_argument("--gzipped", action="store_true", dest="gzipped", default=False, help="input_fasta is gzipped")
parser.add_argument('--output', action='store', dest='output', help='Output file')

args = parser.parse_args()

generate_coverage_table(args.input_fasta, args.input_tabular, args.gzipped, args.output)
