#!/usr/bin/env python

import sys

assert sys.version_info[:2] >= (2, 6)


def __main__():
    skipped_lines = 0
    first_skipped_line = None
    # was sys.argv[2] but we need stdout for a pipe in bam_bed_gff_to_bigwig.xml
    for i, line in enumerate(sys.stdin):
        line = line.rstrip("\r\n")
        if line and not line.startswith("#"):
            try:
                elems = line.split("\t")
                start = str(int(elems[3]) - 1)
                endoff = str(int(elems[4]) - 1)
                # GFF format: chrom, source, name, chromStart, chromEnd, score, strand
                # bedtools puts out only 4 fields: chrom, chromStart, chromEnd, score
                sys.stdout.write(f"{elems[0]}\t{start}\t{endoff}\t0\n")
            except Exception:
                skipped_lines += 1
                if not first_skipped_line:
                    first_skipped_line = i + 1
        else:
            skipped_lines += 1
            if not first_skipped_line:
                first_skipped_line = i + 1


if __name__ == "__main__":
    __main__()
