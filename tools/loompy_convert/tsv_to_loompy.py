import numpy as np
import loompy
import argparse
import os
parser = argparse.ArgumentParser(description="Import tsv files")
parser.add_argument('-f', '--files', help="Tsv file to import")
parser.add_argument('-r', '--rownames', action='store_true', help="Flag for if tsv file contains row names")
parser.add_argument('-c', '--colnames', action='store_true', help="Flag for if tsv file contains column names")
parser.add_argument('-o', '--output', help="Output file path")
parser.add_argument('--row_attr', help="Specifies the title of the row attributes for loompy")
parser.add_argument('--col_attr', help="Specifies the title of the column attributes for loompy")
args = parser.parse_args()
if args.colnames:
    colnames = True
else:
    colnames = False
if args.rownames:
    rownames = True
else:
    rownames = False
if args.row_attr:
    row_info = args.row_attr
else:
    row_info = "Gene"
if args.col_attr:
    col_info = args.col_attr
else:
    col_info = "Cell"
data = args.files
outfile = args.output
if os.path.isfile(outfile) == True:
    os.remove(outfile)
rows = []
matrix = []
with open(data) as d:
    if colnames == True:  #Generate column names from first row
        cols = d.readline()
        cols = cols.strip().split()
    i = 0
    for line in d: #Start ading rows to matrix
        line = list(line.strip().split())
        linelen = len(line)
        if rownames == True: #First element in rows becomes row name
            rows.append(line[0])
            for x in line[1:]:
                matrix.append(float(x))
        else:
            for x in line: #No row name specified, row names omitted
                matrix.append(float(x))
            i += 1
            rows = range(0, i) #Generates row names as numerical sequence from 0 - number of rows - 1
    if colnames == False: #Generates row names as numerical sequence from 0 - number of columns - 1
        cols = list(range(0,int(len(matrix)/len(rows))))
    if colnames == True and rownames == True:
        if len(cols) == linelen:
            raise Exception("Number of column labels incorrect for number of columns. Number of column labels must be one fewer than the total number of columns, as row names are the first column.")
matrix = np.array(matrix).reshape(len(rows), len(cols)) #convert main matrix into numpy array
row_attrs = {row_info: np.array(rows)} #Creates row attribute dictionary for loompy
col_attrs = {col_info: np.array(cols)} #Creates column attribute dictionary for loompy
loompy.create(outfile, matrix, row_attrs, col_attrs) #Creates loompy object from matrix, row attributes and column attributes. Outputs them as outfile.
