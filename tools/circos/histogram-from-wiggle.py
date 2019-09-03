#!/usr/bin/env python
import logging
import sys

import pyBigWig

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

bins = {}

if __name__ == "__main__":
    files = sys.argv[1:]

    for idx, item in enumerate(files):
        bw = pyBigWig.open(item)
        for chrom in bw.chroms().keys():
            for (start, end, value) in bw.intervals(chrom):
                key = "%s\t%s\t%s" % (chrom, start, end)

                if key not in bins:
                    bins[key] = [0] * len(files)

                bins[key][idx] = value

    for (k, v) in bins.items():
        sys.stdout.write(k)
        sys.stdout.write("\t")
        sys.stdout.write(",".join(map(str, v)))
        sys.stdout.write("\n")
