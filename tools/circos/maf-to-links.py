from Bio import AlignIO
import itertools
import sys


for aln in AlignIO.parse(sys.argv[1], "maf"):

    for (a, b) in itertools.combinations(aln, 2):
        a_s = a.annotations['start']
        a_l = a.annotations['size']
        b_s = b.annotations['start']
        b_l = b.annotations['size']

        sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (a.id, a_s, a_s + a_l, b.id, b_s, b_s + b_l))
