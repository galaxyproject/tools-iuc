import sys

import pyBigWig
from Bio import SeqIO
from Bio import SeqUtils


span = int(sys.argv[2])
bw = pyBigWig.open(sys.argv[3], "w")

# Prepare header separately because ugh
data = []
for rec in SeqIO.parse(sys.argv[1], "fasta"):
    data.append((rec.id, len(rec)))
bw.addHeader(data)

for rec in SeqIO.parse(sys.argv[1], "fasta"):
    gc = SeqUtils.GC_skew(rec.seq, span)

    bw.addEntries(rec.id, 0, values=list(gc), span=span, step=span)

bw.close()
