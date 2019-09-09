import itertools
import sys

from Bio import AlignIO


format_mapping = {
    "xmfa": "mauve",
    "maf": "maf",
    "nex": "nexus",
    # 'phylip': 'phylip-relaxed',
    "stockholm": "stockholm",
}

for aln in AlignIO.parse(sys.argv[1], format_mapping.get(sys.argv[2], "maf")):

    for (a, b) in itertools.combinations(aln, 2):
        a_s = a.annotations["start"]
        b_s = b.annotations["start"]

        if "size" in a.annotations:
            a_l = a.annotations["size"]
            b_l = b.annotations["size"]
            a_e = a_l + a_s
            b_e = b_l + b_s
        else:
            a_e = a.annotations["end"]
            b_e = b.annotations["end"]

        sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (a.id, a_s, a_e, b.id, b_s, b_e))
