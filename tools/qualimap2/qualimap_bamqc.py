#!/usr/bin/env python
from __future__ import print_function
import argparse
from subprocess import check_call, CalledProcessError
import shlex
import sys
import logging

log = logging.getLogger(__name__)


def qualimap_bamqc(bam_filename, genomecov_file, out_dir, jv_mem_size):
    cmdline_str = "qualimap bamqc -bam {} -oc {} -outdir {} --java-mem-size={}".format(bam_filename,
                                                                                       genomecov_file,
                                                                                       out_dir,
                                                                                       jv_mem_size)
    cmdline = new_split(cmdline_str)
    try:
        check_call(cmdline)
    except CalledProcessError:
        print("Error running the qualimap bamqc", file=sys.stderr)


def new_split(value):
    lex = shlex.shlex(value)
    lex.quotes = '"'
    lex.whitespace_split = True
    lex.commenters = ''
    return list(lex)


def main():
    parser = argparse.ArgumentParser(description="Generate Bam Quality Statistics")
    parser.add_argument('--input_file')
    parser.add_argument('--out_genome_file', default="bam_stats.genomecov")
    parser.add_argument('--out_dir', default="/tmp/bamstats")
    parser.add_argument('--java_mem_size', default="8G")

    args = parser.parse_args()

    qualimap_bamqc(args.input_file, args.out_genome_file, args.out_dir, args.java_mem_size)


if __name__ == "__main__":
    main()
