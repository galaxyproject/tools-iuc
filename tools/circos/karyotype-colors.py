#!/usr/bin/env python
import sys


highest = int(sys.argv[1])
for idx in range(highest + 1):
    sys.stdout.write(
        "chr{idx}color = lch(66,104,{pos})\n".format(
            idx=idx, pos=int(float(idx) / highest * 360)
        )
    )
