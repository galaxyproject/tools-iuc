#!/usr/bin/env python

import sys


with open(sys.argv[1]) as i:
    bed_data = i.readlines()

sanitized_data = []
try:
    for record in bed_data:
        if record.strip():
            fields = record.split('\t')
            sanitized_data.append(
                '\t'.join(fields[:4] + ['60'] + fields[5:])
            )
except IndexError:
    pass  # leave column number issue to getmasked
else:
    with open(sys.argv[1], 'w') as o:
        o.writelines(sanitized_data)
