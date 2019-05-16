"""This module converts a tsv file into a binary loom file"""

import argparse
import os

import loompy
import numpy as np

parser = argparse.ArgumentParser(description="Loompy file converter flags")
parser.add_argument('--VERSION', action='version', version='%(prog)s 0.1.0',
                    help="Displays tool version")
parser.add_argument('--rowfile', '-r', help="File of row attributes & values")
parser.add_argument('--colfile', '-c',
                    help="File of column attributes and values")
parser.add_argument('--output', '-o', help="Output file name")
parser.add_argument('--files', '-f', nargs='*',
                    help="Input tsv files. First file becomes main layer.")
args = parser.parse_args()

colsfile = args.colfile
rowsfile = args.rowfile
filename = args.output
alldata = args.files

alayers = []
layernames = []
rowdict = {}
coldict = {}

#  Creates dictionary based on row file
#  For each attribute:
#  Attribute: [attribute values]
with open(rowsfile, "r") as rows:
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
#  Same as above, but for columns
with open(colsfile, "r") as cols:
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
#  Finding dimensions for the loom layers
rowshape = len(rowdict[list(rowdict.keys())[0]])
colshape = len(coldict[list(coldict.keys())[0]])

#  Creates a list with each element being entire matrix of
#  each layer file as floats
for file in range(0, len(alldata)):
    layer = alldata[file][:-4]
    layer = layer.split("/")[-1]
    if layer == "":
        raise Exception("Please only use named files")
    layernames.append(layer)
    cfile = alldata[file]
    with open(cfile, "r") as tsv:
        cmatrix = []
        for line in tsv:
            line = line.strip().split("\t")
            line = [float(i) for i in line]
            cmatrix += line
        alayers.append(cmatrix)

#  Loompy cannot overwright existing files. If somehow it finds
#  a second file with the same name, it must be deleted
if os.path.isfile(filename):
    os.remove(filename)
#  To create the file properly, the first row and column attributes must be
#  added separately in the form of individual dictionaries
row_attrs = {row_attributes[0]: np.asarray(rowdict[row_attributes[0]])}
col_attrs = {col_attributes[0]: np.asarray(coldict[col_attributes[0]])}
matrix = np.asarray(alayers[0])
matrix = matrix.astype(float)
matrix = matrix.reshape(rowshape, colshape)
#  Creation of initial loom file
loompy.create(filename, matrix, row_attrs, col_attrs)
#  Adding all row and column attributes, then all layers
with loompy.connect(filename) as loomfile:
    for x in row_attributes:
        loomfile.ra[x] = rowdict[x]
    for y in col_attributes:
        loomfile.ca[y] = coldict[y]
    for z in range(1, len(alayers)):
        matrix = np.asarray(alayers[z])
        matrix = matrix.astype(float)
        matrix = matrix.reshape(rowshape, colshape)
        loomfile[layernames[z]] = matrix
