"""Normalize a BED / region list to a 4-column tab BED (chr, start, end, name).

Reads one file (a user BED, or a Galaxy config file holding typed regions), splits
each non-empty, non-comment line on whitespace, and writes chr, start, end and a
name column (the 4th field if present, otherwise the chromosome). Replaces the awk
one-liners so the tool needs no gawk and never interpolates user text into the shell.
"""
import sys

with open(sys.argv[1]) as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 3:
            continue
        name = f[3] if len(f) >= 4 else f[0]
        sys.stdout.write("\t".join([f[0], f[1], f[2], name]) + "\n")
