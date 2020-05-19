#!/usr/bin/env python
import sys
from Bio import SeqIO


highest = 0

# Process fasta data, extracting only headers
for idx, seq in enumerate(SeqIO.parse(sys.argv[1], "fasta")):
    sys.stdout.write(
        "chr	-	{seq_id}	{seq_id}	0	{length}	chr{idx}color\n".format(
            seq_id=seq.id, idx=idx, length=len(seq)
        )
    )
    highest = idx


with open(sys.argv[2], "w") as handle:
    for idx in range(highest + 1):
        handle.write(
            "chr{idx}color = lch(50,121,{pos})\n".format(
                idx=idx, pos=int(float(idx) / highest * 360)
            )
        )
