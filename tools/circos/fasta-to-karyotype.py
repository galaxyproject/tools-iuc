#!/usr/bin/env python
import sys

from Bio import SeqIO

# Process fasta data, extracting only headers
for idx, seq in enumerate(SeqIO.parse(sys.argv[1], 'fasta')):
    sys.stdout.write("chr - {seq_id} {seq_id} 0 {length} set3-12-qual-{color}\n".format(
        seq_id=seq.id, idx=idx, length=len(seq), color=((idx + 1) % 12)
    ))
