import argparse

import pytrf  # 1.3.0
from pyfastx import Fastx  # 0.5.2

"""
Allows all STR or those for a subset of motifs to be written to a bed file
Designed to build some of the microsatellite tracks from https://github.com/arangrhie/T2T-Polish/tree/master/pattern for the VGP.
"""


def write_ssrs(args):
    """
    The integers in the call change the minimum repeats for mono-, di-, tri-, tetra-, penta-, hexa-nucleotide repeats
    ssrs = pytrf.STRFinder(name, seq, 10, 6, 4, 3, 3, 3)
    NOTE: Dinucleotides GA and AG are reported separately by https://github.com/marbl/seqrequester.
    The reversed pair STRs are about as common in the documentation sample.
    Sequence read bias might be influenced by GC density or some other specific motif.
    """
    bed = []
    specific = None
    if args.specific:
        specific = args.specific.upper().split(",")
    fa = Fastx(args.fasta, uppercase=True)
    for name, seq in fa:
        for ssr in pytrf.STRFinder(
            name,
            seq,
            args.monomin,
            args.dimin,
            args.trimin,
            args.tetramin,
            args.pentamin,
            args.hexamin,
        ):
            row = (
                ssr.chrom,
                ssr.start - 1,
                ssr.end,
                ssr.motif,
                ssr.repeat,
                ssr.length,
            )
            # pytrf reports a 1 based start position so start-1 fixes the bed interval lengths
            if args.specific and ssr.motif in specific:
                bed.append(row)
            elif args.mono and len(ssr.motif) == 1:
                bed.append(row)
            elif args.di and len(ssr.motif) == 2:
                bed.append(row)
            elif args.tri and len(ssr.motif) == 3:
                bed.append(row)
            elif args.tetra and len(ssr.motif) == 4:
                bed.append(row)
            elif args.penta and len(ssr.motif) == 5:
                bed.append(row)
            elif args.hexa and len(ssr.motif) == 6:
                bed.append(row)
    bedtosort = [(x[0], x[1], x[2], x) for x in bed]
    bedtosort.sort()
    # decorate and undecorate to avoid bogus alphanumeric ordering
    sbed = [x[3] for x in bedtosort]
    obed = ["%s\t%d\t%d\t%s_%d\t%d" % x for x in sbed]
    with open(args.bed, "w") as outbed:
        outbed.write("\n".join(obed))
        outbed.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    a = parser.add_argument
    a("--di", action="store_true")
    a("--tri", action="store_true")
    a("--tetra", action="store_true")
    a("--penta", action="store_true")
    a("--hexa", action="store_true")
    a("--mono", action="store_true")
    a("--dimin", default=2, type=int)
    a("--trimin", default=2, type=int)
    a("--tetramin", default=2, type=int)
    a("--pentamin", default=2, type=int)
    a("--hexamin", default=2, type=int)
    a("--monomin", default=2, type=int)
    a("-f", "--fasta", default="humsamp.fa")
    a("-b", "--bed", default="humsamp.bed")
    a("--specific", default=None)
    a("--minreps", default=2, type=int)
    args = parser.parse_args()
    write_ssrs(args)
