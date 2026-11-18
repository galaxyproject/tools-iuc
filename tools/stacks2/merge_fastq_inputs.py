#!/usr/bin/env python3
"""Merge FASTQ/FASTA input files into a single file.

Since Stacks 2.65 process_shortreads can no longer scan a directory for input
files (only process_radtags enables that mode), hence all input files have to
be passed via -f (or -1/-2 for paired-end data).  This script concatenates an
arbitrary number of (optionally gzipped) input files into one file.

usage: merge_fastq_inputs.py OUTPUT FILE [FILE ...]

The record size (4 lines for FASTQ, 2 lines for FASTA) is determined from the
first input file.
"""
import gzip
import sys


def open_input(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def open_output(path):
    if path.endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "w")


def record_size(path):
    with open_input(path) as fh:
        first = fh.read(1)
    return 2 if first == ">" else 4


def copy_records(path, out, size):
    with open_input(path) as fh:
        while True:
            lines = [fh.readline() for _ in range(size)]
            if not lines[0]:
                return
            out.write("".join(line if line.endswith("\n") else line + "\n" for line in lines))


def main(argv):
    if len(argv) < 3:
        sys.stderr.write(__doc__)
        return 1
    out_path, infiles = argv[1], argv[2:]
    size = record_size(infiles[0])
    with open_output(out_path) as out:
        for path in infiles:
            copy_records(path, out, size)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
