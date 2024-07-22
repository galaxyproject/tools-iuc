import argparse
import subprocess

import pytrf  # 1.3.0
from pyfastx import Fastx  # 0.5.2

"""
Allows all STR or those for a subset of motifs to be written to a bed file
Designed to build some of the microsatellite tracks from https://github.com/arangrhie/T2T-Polish/tree/master/pattern for the VGP.
"""


def getDensity(name, bed, chrlen, winwidth):
    """
    pybigtools can write bigwigs and they are processed by other ucsc tools - but jb2 will not read them.
    Swapped the conversion to use a bedgraph file processed by bedGraphToBigWig
    """
    nwin = int(chrlen / winwidth)
    d = [0.0 for x in range(nwin + 1)]
    for b in bed:
        nt = b[5]
        bin = int(b[1] / winwidth)
        d[bin] += nt
    bedg = [
        (name, (x * winwidth), ((x + 1) * winwidth) - 1, float(d[x]))
        for x in range(nwin + 1)
        if (x + 1) * winwidth <= chrlen
    ]
    return bedg


def write_ssrs(args):
    """
    The integers in the call change the minimum repeats for mono-, di-, tri-, tetra-, penta-, hexa-nucleotide repeats
    ssrs = pytrf.STRFinder(name, seq, 10, 6, 4, 3, 3, 3)
    NOTE: Dinucleotides GA and AG are reported separately by https://github.com/marbl/seqrequester.
    The reversed pair STRs are about as common in the documentation sample.
    Sequence read bias might be influenced by GC density or some other specific motif.
    """
    bed = []
    wig = []
    chrlens = {}
    specific = None
    if args.specific:
        specific = args.specific.upper().split(",")
    fa = Fastx(args.fasta, uppercase=True)
    for name, seq in fa:
        chrlen = len(seq)
        chrlens[name] = chrlen
        cbed = []
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
                ssr.start,
                ssr.end,
                ssr.motif,
                ssr.repeat,
                ssr.length,
            )
            if args.specific and ssr.motif in specific:
                cbed.append(row)
            elif args.mono and len(ssr.motif) == 1:
                cbed.append(row)
            elif args.di and len(ssr.motif) == 2:
                cbed.append(row)
            elif args.tri and len(ssr.motif) == 3:
                cbed.append(row)
            elif args.tetra and len(ssr.motif) == 4:
                cbed.append(row)
            elif args.penta and len(ssr.motif) == 5:
                cbed.append(row)
            elif args.hexa and len(ssr.motif) == 6:
                cbed.append(row)
        if args.bigwig:
            w = getDensity(name, cbed, chrlen, args.winwidth)
            wig += w
        bed += cbed
    if args.bigwig:
        wig.sort()
        bedg = ["%s %d %d %.2f" % x for x in wig]
        with open("temp.bedg", "w") as bw:
            bw.write("\n".join(bedg))
        chroms = ["%s\t%s" % (x, chrlens[x]) for x in chrlens.keys()]
        with open("temp.chromlen", "w") as cl:
            cl.write("\n".join(chroms))
        cmd = ["bedGraphToBigWig", "temp.bedg", "temp.chromlen", args.bed]
        subprocess.run(cmd)
    else:
        bed.sort()
        obed = ["%s\t%d\t%d\t%s_%d\t%d" % x for x in bed]
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
    a("--bigwig", action="store_true")
    a("--winwidth", default=128, type=int)
    a("--specific", default=None)
    a("--minreps", default=2, type=int)
    args = parser.parse_args()
    write_ssrs(args)
