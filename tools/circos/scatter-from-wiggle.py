#!/usr/bin/env python
import pyBigWig
import logging
import sys

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

if __name__ == '__main__':
    bw = pyBigWig.open(sys.argv[1])
    for chrom in bw.chroms().keys():
        for (start, end, value) in bw.intervals(chrom):
            sys.stdout.write('%s\t%s\t%s\t%s\n' % (chrom, start, end, value))
