#!/usr/bin/env python

###################################################################
#
# gbk2orf.py by Errol Strain (estrain@gmail.com)
#
# Read a GenBank file and export fasta formatted amino acid and
# CDS files
#
###################################################################

import sys
from optparse import OptionParser

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Command line usage
usage = "usage: %prog -g input.gbk -a aa.fasta -n nuc.fasta"
p = OptionParser(usage)
p.add_option("-t", "--translate", dest="transtabl", type="int", default=11,
             help="Translation table used to translate coding regions (default=11)")
p.add_option("-g", "--genbank", dest="gb_file", help="GenBank input file")
p.add_option("-a", "--amino_acid", dest="aa_file", help="Fasta amino acid output")
p.add_option("-n", "--nucleotide", dest="orf_file", help="Fasta nucleotide output")
(opts, args) = p.parse_args()
# Do I need this next line?
if not opts and not args:
    p.error("Use --help to see usage")
if len(sys.argv) == 1:
    p.error("Use --help to see usage")

# Lists to hold SeqRecords
aalist = []
nuclist = []

# If the CDS does not have a locus tag the name will be assigned using the
# order in which it was found
feat_count = 0

# Iterate through genbank records in input file
for gb_record in SeqIO.parse(open(opts.gb_file, "r"), "genbank"):
    for (index, feature) in enumerate(gb_record.features):
        if feature.type == "CDS":
            feat_count = feat_count + 1
            gene = feature.extract(gb_record.seq)
            if "locus_tag" in feature.qualifiers:
                value = feature.qualifiers["locus_tag"][0]
            else:
                value = "Index_" + str(feat_count)
            nuclist.append(SeqRecord(Seq(str(gene)), id=value, name=value))
            pro = Seq(str(gene.translate(table=opts.transtabl, to_stop=True)))
            aalist.append(SeqRecord(pro, id=value, name=value))

# Write out lists in fasta format
aa_handle = open(opts.aa_file, "w")
SeqIO.write(aalist, aa_handle, "fasta")
aa_handle.close()
orf_handle = open(opts.orf_file, "w")
SeqIO.write(nuclist, orf_handle, "fasta")
orf_handle.close()
