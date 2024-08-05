#!/usr/bin/env python
import logging
import sys

import pyBigWig

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

if __name__ == "__main__":
    files = sys.argv[1:]
    bws = [pyBigWig.open(x) for x in files]

    # obtain some chroms. Hope all sets are identical!
    k = list(bws[0].chroms().keys())

    # do magic?
    # nah.
    # just assert that intervals are identical.
    # and crash otherwise.
    # sorry not sorry.

    for chrom in k:
        for interval_set in zip(*[bw.intervals(chrom) for bw in bws]):
            (start, end) = interval_set[0][0:2]
            values = ",".join(map(str, [x[2] for x in interval_set]))
            sys.stdout.write("%s\t%s\t%s\t%s\n" % (chrom, start, end, values))
