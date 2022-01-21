#!/usr/bin/env python

"""Converts a loompy file to tsv file(s). Each layer becomes a new file."""

import argparse

import loompy

parser = argparse.ArgumentParser(description="Loompy file converter flags")
parser.add_argument('--version', action='version', version='%(prog)s 0.1.0',
                    help="Displays tool version")
parser.add_argument("-f", "--file", help="loom file to import")
args = parser.parse_args()

file = args.file

matrices = []
allcols = []
colstrings = []
allrows = []

# Build background info for all attributes and layers
loompyfile = loompy.connect(file)
row_attributes = loompyfile.ra.keys()  # List of row attributes
for row in row_attributes:  # Each list represents rownames for row_attributes
    c_row = loompyfile.ra[row]
    c_row = [str(r) for r in c_row]
    allrows.append(c_row)
col_attributes = loompyfile.ca.keys()  # List of column attributes
for col in col_attributes:  # each list represents colnames for col_attributes
    c_col = loompyfile.ca[col]
    c_col = [str(c) for c in c_col]
    allcols.append(c_col)
layers = loompyfile.layers.keys()  # List of layers
for layer in layers:  # List with each element being a loompy layer
    c_layer = loompyfile[layer]
    c_layer = c_layer[:, :]
    c_layer = c_layer.astype(str)
    matrices.append(c_layer)

# Create column attribute output
with open("attributes/col_attr.tsv", "w") as colout:
    col_attributes = "\t".join(col_attributes) + "\n"
    colout.write(col_attributes)
    for length in range(0, len(c_col)):
        attributestring = ""
        for col in allcols:
            attributestring = attributestring + col[length] + "\t"
        while attributestring[-1] == "\t":
            attributestring = attributestring[:-1]
        colout.write(attributestring)
        colout.write("\n")
# Create row attribute output
with open("attributes/row_attr.tsv", "w") as rowout:
    row_attributes = "\t".join(row_attributes) + "\n"
    rowout.write(row_attributes)
    for length in range(0, len(c_row)):
        attributestring = ""
        for row in allrows:
            attributestring = attributestring + row[length] + "\t"
        while attributestring[-1] == "\t":
            attributestring = attributestring[:-1]
        rowout.write(attributestring)
        rowout.write("\n")

# Build output files for each layer
for x in range(0, len(layers)):
    # Output file name generation
    if layers[x] in layers[0: x]:  # Different output names if layers have same names somehow
        repeats = layers[0, x].count(layer[x])
        outputname = "output/" + layers[x] + repeats + ".tsv"
    elif layers[x] == "":  # Empty layer name
        outputname = "output/mainmatrix.tsv"
    else:
        outputname = "output/" + str(layers[x]) + ".tsv"  # Usual case
# Matrix output
    with open(outputname, "w") as outputmatrix:
        for line in matrices[x]:
            line = "\t".join(line)
            line += "\n"
            line = line
            outputmatrix.write(line)
