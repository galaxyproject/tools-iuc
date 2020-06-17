#!/usr/bin/env python
"""This program adds layers, row attributes or column attributes for loom files"""

import argparse

import loompy
import numpy as np

parser = argparse.ArgumentParser(description="Loompy file converter flags")
parser.add_argument('--VERSION', action='version', version='%(prog)s 0.1.0',
                    help="Displays tool version")
parser.add_argument('--file', '-f',
                    help="Loom file to which data will be added")
parser.add_argument('--rowfile', '-r', help="File of row attributes & values")
parser.add_argument('--colfile', '-c',
                    help="File of column attributes and values")
parser.add_argument('--layers', '-l', nargs='*',
                    help="Input tsv files. First file becomes main layer.")
parser.add_argument('--add', '-a', choices=["rows", "cols", "layers"],
                    help="Selects rows, columns or layers to be added to file")
args = parser.parse_args()

lfile = args.file
if args.rowfile:
    rowfile = args.rowfile
if args.colfile:
    colfile = args.colfile
if args.layers:
    alllayers = args.layers
addselect = args.add
# Check proper flags for chosen attributes are being added
if addselect == "cols" and not args.colfile:
    raise Exception("To add column attributes, column flag and file must be provided")
if addselect == "rows" and not args.rowfile:
    raise Exception("To add row attributes, row flag and file must be provided")
if addselect == "layers" and not args.layers:
    raise Exception("To add layers, a layer flag and file(s) must be provided")

layernames = []
rowdict = {}
coldict = {}

with loompy.connect(lfile) as loomfile:
    # Loom file dimensions
    nrow = loomfile.shape[0]
    ncol = loomfile.shape[1]
    if addselect == "layers":
        layernames = []
        # Generate layer names based on file names
        for x in range(0, len(alllayers)):
            layer = alllayers[x]
            layer = layer.split("/")[-1].split(".")[-2]  # Takes away path, takes off extension
            layernames.append(layer)
        # Add in the layers themselves
        for layer in range(0, len(alllayers)):
            matrix = ""
            with open(alllayers[layer], "r") as infile:
                rows = 0
                count = 0
                for line in infile:
                    if count == 0:
                        cols = len(line.split("\t"))
                        if cols != ncol:
                            raise Exception("Dimensions of new matrix incorrect for this loom file. New matrices must be %d by %d" % (nrow, ncol))
                    matrix = matrix + line + "\t"
                    rows += 1
                if rows != nrow:
                    raise Exception("Dimensions of new matrix incorrect for this loom file. New matrices must be %d by %d")
            matrix = matrix.split("\t")
            matrix = [float(n) for n in matrix[:-1]]
            matrix = np.asarray(matrix).reshape(nrow, ncol)
            loomfile[layernames[layer]] = matrix
    elif addselect == "rows":
        with open(rowfile, "r") as rows:
            count = 0
            for line in rows:
                line = line.strip().split("\t")
                if count == 0:  # First time through
                    row_attributes = line
                    for x in row_attributes:
                        rowdict[x] = []
                    count += 1
                else:
                    for x in range(0, len(line)):
                        rowdict[row_attributes[x]].append(line[x])
        for x in row_attributes:
            if len(rowdict[x]) != nrow:
                raise Exception("Incorrect length of row. Row length must be: %d" % nrow)
            loomfile.ra[x] = rowdict[x]
    elif addselect == "cols":
        with open(colfile, "r") as cols:
            count = 0
            for line in cols:
                line = line.replace('\"', "")
                line = line.replace(' ', "")
                line = line.strip().split("\t")
                if count == 0:  # First time through
                    col_attributes = line
                    for x in col_attributes:
                        coldict[x] = []
                    count += 1
                else:
                    for x in range(0, len(line)):
                        coldict[col_attributes[x]].append(line[x])
        for y in col_attributes:
            if len(coldict[y]) != ncol:
                raise Exception("Incorrect length of column. Column length must be: %d" % ncol)
            loomfile.ca[y] = coldict[y]
