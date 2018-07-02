#!/usr/bin/env python
import wiggle
import logging
import sys

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

if __name__ == '__main__':
    walker = wiggle.Wiggle()
    for chrom, start, end, value in walker.walk(sys.stdin):
        print('%s %s %s %s' % (chrom, start, end, value))
