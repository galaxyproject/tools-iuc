"""Converts a loompy file to tsv file(s). Each layer becomes a new file."""

import argparse

import loompy

PARSER = argparse.ArgumentParser(description="Loompy file converter flags")
PARSER.add_argument('--VERSION', action='version', version='%(prog)s 0.1.0',
                    help="Displays tool version")
PARSER.add_argument("-f", "--FILE", help="loom file to import")
ARGS = PARSER.parse_args()

FILE = ARGS.FILE

MATRICES = []
ALLCOLS = []
COLSTRINGS = []
ALLROWS = []

# Build background info for all attributes and layers
LOOMPYFILE = loompy.connect(FILE)
ROW_ATTRIBUTES = LOOMPYFILE.ra.keys()  # List of row attributes
for ROW in ROW_ATTRIBUTES:  # Each list represents rownames for row_attributes
    C_ROW = LOOMPYFILE.ra[ROW]
    C_ROW = [str(R) for R in C_ROW]
    ALLROWS.append(C_ROW)
COL_ATTRIBUTES = LOOMPYFILE.ca.keys()  # List of column attributes
for COL in COL_ATTRIBUTES:  # Each list represents colnames for col_attributes
    C_COL = LOOMPYFILE.ca[COL]
    C_COL = [str(C) for C in C_COL]
    ALLCOLS.append(C_COL)
LAYERS = LOOMPYFILE.layers.keys()  # List of layers
for LAYER in LAYERS:  # List with each element being a loompy layer
    C_LAYER = LOOMPYFILE[LAYER]
    C_LAYER = C_LAYER[:, :]
    C_LAYER = C_LAYER.astype(str)
    MATRICES.append(C_LAYER)

# Create column attribute output
with open("attributes/col_attr.tsv", "w") as COLOUT:
    COL_ATTRIBUTES = "\t".join(COL_ATTRIBUTES) + "\n"
    COLOUT.write(COL_ATTRIBUTES)
    for LENGTH in range(0, len(C_COL)):
        ATTRIBUTESTRING = ""
        for COL in ALLCOLS:
            ATTRIBUTESTRING = ATTRIBUTESTRING + COL[LENGTH] + "\t"
        ATTRIBUTESTRING = ATTRIBUTESTRING[:-1]
        COLOUT.write(ATTRIBUTESTRING)
        COLOUT.write("\n")
# Create row attribute output
with open("attributes/row_attr.tsv", "w") as ROWOUT:
    ROW_ATTRIBUTES = "\t".join(ROW_ATTRIBUTES) + "\n"
    ROWOUT.write(ROW_ATTRIBUTES)
    for LENGTH in range(0, len(C_ROW)):
        ATTRIBUTESTRING = ""
        for ROW in ALLROWS:
            ATTRIBUTESTRING = ATTRIBUTESTRING + ROW[LENGTH] + "\1"
        ATTRIBUTESTRING = ATTRIBUTESTRING[:-1]
        ROWOUT.write(ATTRIBUTESTRING)
        ROWOUT.write("\n")

# Build output files for each layer
for X in range(0, len(LAYERS)):
    # Output file name generation
    if LAYERS[X] in LAYERS[0: X]:  # Different output names if layers have same names somehow
        REPEATS = LAYERS[0, X].count(LAYER[X])
        OUTPUTNAME = "output/" + LAYERS[X] + REPEATS + ".tsv"
    elif LAYERS[X] == "":  # Empty layer name
        OUTPUTNAME = "output/mainmatrix.tsv"
    else:
        OUTPUTNAME = "output/" + str(LAYERS[X]) + ".tsv"  # Usual case
# Matrix output
    with open(OUTPUTNAME, "w") as OUTPUTMATRIX:
        for LINE in MATRICES[X]:
            LINE = "\t".join(LINE)
            LINE += "\n"
            OUTPUTMATRIX.write(LINE)
