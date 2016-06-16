#!/usr/bin/env python
from Bio import SeqIO
import sys

for idx, seq in enumerate(SeqIO.parse(sys.argv[1], 'fasta')):
    print "chr - {seq_id} {idx} 0 {length} set3-12-qual-{color}\n".format(
        seq_id=seq.id, idx=idx, length=len(seq), color=((idx + 1) % 12)
    )
