#!/usr/bin/env python
import logging
import wiggle
import sys

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

bins = {}

if __name__ == '__main__':
    walker = wiggle.Wiggle()
    files = sys.argv[1:]

    for idx, item in enumerate(files):
        with open(item, 'r') as handle:
            for chrom, start, end, value in walker.walk(handle):
                key = '%s %s %s' % (chrom, start, end)
                if key not in bins:
                    bins[key] = [0] * len(files)
                bins[key][idx] = value

    for (k, v) in bins.items():
        sys.stdout.write(k)
        sys.stdout.write(' ')
        sys.stdout.write(' '.join(map(str, v)))
        sys.stdout.write('\n')
