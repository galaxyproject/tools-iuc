"""
Ross Lazarus June 2024 for VGP
Bigwigs are great, but hard to reliably "see" small low coverage or small very high coverage regions.
Colouring in JB2 tracks will need a new plugin, so this code will find bigwig regions above and below a chosen percentile point.
0.99 and 0.01 work well in testing with a minimum span of 10 bp.
Multiple bigwigs **with the same reference** can be combined - bed segments will be named appropriately
Combining multiple references works but is silly because only display will rely on one reference so others will not be shown...
Tricksy numpy method from http://gregoryzynda.com/python/numpy/contiguous/interval/2019/11/29/contiguous-regions.html
takes about 95 seconds for a 17MB test wiggle
JBrowse2 bed normally displays ignore the score, so could provide separate low/high bed file outputs as an option.
Update june 30 2024: wrote a 'no-build' plugin for beds to display red/blue if >0/<0 so those are used for scores
Bed interval naming must be short for JB2 but needs input bigwig name and (lo or hi).
"""

import argparse
import copy
import os
import sys
from pathlib import Path

import numpy as np
import pybigtools


class findOut:

    def __init__(self, args):
        self.bwnames = args.bigwig
        self.bwlabels = args.bigwiglabels
        self.bedwin = args.minwin
        self.outbeds = args.outbeds
        self.bedouthi = args.bedouthi
        self.bedoutlo = args.bedoutlo
        self.bedouthilo = args.bedouthilo
        self.tableoutfile = args.tableoutfile
        self.bedwin = args.minwin
        self.qhi = args.qhi
        self.qlo = None
        try:
            f = float(args.qlo)
            self.qlo = f
        except:
            print('qlo=',args.qlo)
        nbw = len(args.bigwig)
        nlab = len(args.bigwiglabels)
        if nlab < nbw:
            self.bwlabels += ["Nolabel"] * (nbw - nlab)
        self.makeBed()

    def processVals(self, bw, isTop):
        """
        idea from http://gregoryzynda.com/python/numpy/contiguous/interval/2019/11/29/contiguous-regions.html
        Fast segmentation into regions by taking np.diff on the boolean array of over (under) cutpoint indicators in bwex.
        This only gives non-zero values at the segment boundaries where there's a change, so those zeros are all removed in bwexdnz
        leaving an array of segment start/end positions. That's twisted around into an array of start/end coordinates.
        Magical. Fast. Could do the same for means or medians over windows for sparse bigwigs like repeat regions.
        """
        if isTop:
            bwex = np.r_[False, bw >= self.bwtop, False]  # extend with 0s
        else:
            bwex = np.r_[False, bw <= self.bwbot, False]
        bwexd = np.diff(bwex)
        bwexdnz = bwexd.nonzero()[0]
        bwregions = np.reshape(bwexdnz, (-1, 2))
        return bwregions

    def writeBed(self, bed, bedfname):
        """
        potentially multiple
        """
        bed.sort()
        beds = ["%s\t%d\t%d\t%s\t%d" % x for x in bed]
        with open(bedfname, "w") as bedf:
            bedf.write("\n".join(beds))
            bedf.write("\n")

    def makeTableRow(self, bw, bwlabel, chr):
        """
        called for every contig, but messy inline
        """
        bwmean = np.mean(bw)
        bwstd = np.std(bw)
        bwmax = np.max(bw)
        nrow = np.size(bw)
        bwmin = np.min(bw)
        row = "%s\t%s\t%d\t%f\t%f\t%f\t%f" % (
            bwlabel,
            chr,
            nrow,
            bwmean,
            bwstd,
            bwmin,
            bwmax,
        )
        if self.qhi is not None:
            row += "\t%.2f" % self.bwtop
        else:
            row += "\tnoqhi"
        if self.qlo is not None:
            row += "\t%.2f" % self.bwbot
        else:
            row += "\tnoqlo"
        return row

    def makeBed(self):
        bedhi = []
        bedlo = []
        bwlabels = self.bwlabels
        bwnames = self.bwnames
        if self.tableoutfile:
            restab = ["bigwig\tcontig\tn\tmean\tstd\tmin\tmax\tqtop\tqbot"]
        for i, bwname in enumerate(bwnames):
            bwlabel = bwlabels[i].replace(" ", "")
            fakepath = "in%d.bw" % i
            if os.path.isfile(fakepath):
                os.remove(fakepath)
            p = Path(fakepath)
            p.symlink_to(bwname)  # required by pybigtools (!)
            bwf = pybigtools.open(fakepath)
            chrlist = bwf.chroms()
            chrs = list(chrlist.keys())
            chrs.sort()
            for chr in chrs:
                bw = bwf.values(chr)
                bw = bw[~np.isnan(bw)]  # some have NaN if parts of a contig not covered
                if self.qhi is not None:
                    self.bwtop = np.quantile(bw, self.qhi)
                    bwhi = self.processVals(bw, isTop=True)
                    for j, seg in enumerate(bwhi):
                        if seg[1] - seg[0] >= self.bedwin:
                            score = np.sum(bw[seg[0]:seg[1]])
                            bedhi.append(
                                (
                                    chr,
                                    seg[0],
                                    seg[1],
                                    "%s_%d" % (bwlabel, score),
                                    score,
                                )
                            )
                if self.qlo is not None:
                    self.bwbot = np.quantile(bw, self.qlo)
                    bwlo = self.processVals(bw, isTop=False)
                    for j, seg in enumerate(bwlo):
                        if seg[1] - seg[0] >= self.bedwin:
                            score = -1 * np.sum(bw[seg[0]:seg[1]])
                            bedlo.append(
                                (
                                    chr,
                                    seg[0],
                                    seg[1],
                                    "%s_%d" % (bwlabel, score),
                                    score,
                                )
                            )
                if self.tableoutfile:
                    row = self.makeTableRow(bw, bwlabel, chr)
                    restab.append(copy.copy(row))
        if self.tableoutfile:
            stable = "\n".join(restab)
            with open(self.tableoutfile, "w") as t:
                t.write(stable)
                t.write("\n")
        some = False
        if self.qlo:
            if self.outbeds in ["outall", "outlo", "outlohi"]:
                self.writeBed(bedlo, self.bedoutlo)
                some = True
        if self.qhi:
            if self.outbeds in ["outall", "outlohi", "outhi"]:
                self.writeBed(bedhi, self.bedouthi)
                some = True
        if self.outbeds in ["outall", "outhilo"]:
            allbed = bedlo + bedhi
            self.writeBed(allbed, self.bedouthilo)
            some = True
        if not some:
            sys.stderr.write(
                "Invalid configuration - no output could be created. Was qlo missing and only low output requested for example?"
            )
            sys.exit(2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    a = parser.add_argument
    a("-m", "--minwin", default=10, type=int)
    a("-l", "--qlo", default=None)
    a("-i", "--qhi", default=None, type=float)
    a("--bedouthi", default=None)
    a("--bedoutlo", default=None)
    a("--bedouthilo", default=None)
    a("-w", "--bigwig", nargs="+")
    a("-n", "--bigwiglabels", nargs="+")
    a("-o", "--outbeds", default="outhilo", help="optional high and low combined bed")
    a("-t", "--tableoutfile", default=None)
    args = parser.parse_args()
    findOut(args)
