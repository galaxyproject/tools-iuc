#!/usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('bcfile', help='barcode file')
args = parser.parse_args()

barcodes = []

with open(args.bcfile, "r") as fh:
    for line in fh:
        if len(line) == 0:
            continue
        if line.startswith("#"):
            continue
        barcodes.append(line.split())

if len(barcodes) <= 1:
    sys.exit("barcode file is empty")

# check that all lines have the same number of columns
ncol = None
for bc in barcodes:
    if ncol is None:
        ncol = len(bc)
    elif ncol != len(bc):
        sys.exit("barcode file has inconsistent number of columns")

isname = False
for bc in barcodes:
    if len(bc[-1].strip("ATCGatcg")) > 0:
        isname = True
        break

names = set()
for bc in barcodes:
    if isname:
        n = bc[-1]
    else:
        n = '-'.join(bc)
    if n in names:
        sys.exit("duplicate sample %s in barcode file" % n)
    names.add(n)
