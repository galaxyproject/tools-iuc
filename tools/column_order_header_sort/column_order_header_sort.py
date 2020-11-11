#!/usr/bin/env python

import sys

input_filename = sys.argv[1]
output_filename = sys.argv[2]
delimiter = sys.argv[3]
key_column = sys.argv[4]

try:
    key_column = int(key_column) - 1
except Exception:
    key_column = None

header = None
with open(input_filename, 'r') as fh:
    header = fh.readline().strip('\r\n')
header = header.split(delimiter)
assert len(header) == len(set(header)), "Header values must be unique"
sorted_header = list(header)
if key_column is None:
    columns = []
else:
    columns = [key_column]
    sorted_header.pop(key_column)
sorted_header.sort()

for key in sorted_header:
    columns.append(header.index(key))

with open(input_filename, 'r') as in_fh, open(output_filename, 'w') as out_fh:
    for line in in_fh:
        line = line.strip('\r\n')
        line = line.split(delimiter)

        out_fh.write(delimiter.join([line[j] for j in columns]) + "\n")
