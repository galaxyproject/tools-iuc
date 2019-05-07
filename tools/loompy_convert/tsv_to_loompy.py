"""This module converts a tsv file into a binary loom file"""

import argparse

import os

import loompy

import numpy as np

PARSER = argparse.ArgumentParser(description="Import tsv files")
PARSER.add_argument('-f', '--FILES', help="Tsv file to import")
PARSER.add_argument('-r', '--ROWNAMES', action='store_true', help="tsv file contains row names")
PARSER.add_argument('-c', '--COLNAMES', action='store_true', help="tsv file contains column names")
PARSER.add_argument('-o', '--OUTPUT', help="Output file path")
PARSER.add_argument('--ROW_ATTR', help="Specifies the title of the row attributes for loompy")
PARSER.add_argument('--COL_ATTR', help="Specifies the title of the column attributes for loompy")
ARGS = PARSER.parse_args()
COLNAMES = bool(ARGS.COLNAMES)
ROWNAMES = bool(ARGS.ROWNAMES)
ROW_ATTR = bool(ARGS.ROW_ATTR)
if ARGS.ROW_ATTR:
    ROW_INFO = ARGS.ROW_ATTR
else:
    ROW_INFO = "Gene"
if ARGS.COL_ATTR:
    COL_INFO = ARGS.COL_ATTR
else:
    COL_INFO = "Cell"
DATA = ARGS.FILES
OUTFILE = ARGS.OUTPUT
if os.path.isfile(OUTFILE):
    os.remove(OUTFILE)
ROWS = []
MATRIX = []
with open(DATA) as D:
# Generate column names from first row
    if COLNAMES:
        COLS = D.readline()
        COLS = COLS.strip().split()
    i = 0
    for LINE in D:  # Start ading rows to matrix
        LINE = list(LINE.strip().split())
        LINELEN = len(LINE)
        if ROWNAMES:  # First element in rows becomes row name
            ROWS.append(LINE[0])
            for x in LINE[1:]:
                MATRIX.append(float(x))
        else:
            for x in LINE:  # No row name specified, row names omitted
                MATRIX.append(float(x))
            i += 1
            ROWS = range(0, i)  # Generates row names as numerical sequence from 0 - number of rows minus 1
    if not COLNAMES:
        COLS = list(range(0, int(len(MATRIX)/len(ROWS))))
    if COLNAMES and ROWNAMES:
        if len(COLS) == LINELEN:
            raise Exception("Number of column labels incorrect for number of columns. Number of column labels must be one fewer than the total number of columns, as row names are the first column.")
# Convert main matrix into numpy array
MATRIX = np.array(MATRIX).reshape(len(ROWS), len(COLS))
# Creates row attribute dictionary for loompy
ROW_ATTRS = {ROW_INFO: np.array(ROWS)}
# Creates column attribute dictionary for loompy
COL_ATTRS = {COL_INFO: np.array(COLS)}
# Creates loompy object from matrix, row attributes and column attributes. Outputs them as outfile.
loompy.create(OUTFILE, MATRIX, ROW_ATTRS, COL_ATTRS)
