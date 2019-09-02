#!/usr/bin/env python
import csv
import sys


idx = 0
with open(sys.argv[1], 'r') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
    for row in spamreader:
        if len(row) < 2:
            continue

        seq_id = row[0]
        length = row[1]

        sys.stdout.write("chr - {seq_id} {seq_id} 0 {length} set3-12-qual-{color}\n".format(
            seq_id=seq_id, idx=idx, length=length, color=((idx + 1) % 12)
        ))
        idx += 1
