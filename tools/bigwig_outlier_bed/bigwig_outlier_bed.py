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
import os
import sys
from pathlib import Path

import numpy as np
import pybigtools


class asciihist:

    def __init__(
        self,
        data,
        bins=10,
        minmax=None,
        str_tag="",
        scale_output=80,
        generate_only=True,
    ):
        """
        https://gist.github.com/bgbg/608d9ef4fd75032731651257fe67fc81
        Create an ASCII histogram from an interable of numbers.
        Author: Boris Gorelik boris@gorelik.net. based on  http://econpy.googlecode.com/svn/trunk/pytrix/pytrix.py
        License: MIT
        """
        self.data = data
        self.minmax = minmax
        self.str_tag = str_tag
        self.bins = bins
        self.generate_only = generate_only
        self.scale_output = scale_output
        self.itarray = np.asanyarray(self.data)
        if self.minmax == "auto":
            self.minmax = np.percentile(data, [5, 95])
            if self.minmax[0] == self.minmax[1]:
                # for very ugly distributions
                self.minmax = None
        if self.minmax is not None:
            # discard values that are outside minmax range
            mn = self.minmax[0]
            mx = self.minmax[1]
            self.itarray = self.itarray[self.itarray >= mn]
            self.itarray = self.itarray[self.itarray <= mx]

    def draw(self):
        values, counts = np.unique(self.data, return_counts=True)
        if len(values) <= 20:
            self.bins = len(values)
        ret = []
        if self.itarray.size:
            total = len(self.itarray)
            counts, cutoffs = np.histogram(self.itarray, bins=self.bins)
            cutoffs = cutoffs[1:]
            if self.str_tag:
                self.str_tag = "%s " % self.str_tag
            else:
                self.str_tag = ""
            if self.scale_output is not None:
                scaled_counts = counts.astype(float) / counts.sum() * self.scale_output
            else:
                scaled_counts = counts
            footerbar = "{:s}{:s} |{:s} |".format(
                self.str_tag,
                "-" * 12,
                "-" * 12,
            )
            if self.minmax is not None:
                ret.append(
                    "Trimmed to range (%s - %s)"
                    % (str(self.minmax[0]), str(self.minmax[1]))
                )
            for cutoff, original_count, scaled_count in zip(
                cutoffs, counts, scaled_counts
            ):
                ret.append(
                    "{:s}{:>12.2f} |{:>12,d} | {:s}".format(
                        self.str_tag, cutoff, original_count, "*" * int(scaled_count)
                    )
                )
            ret.append(footerbar)
            ret.append("{:s}{:>12s} |{:>12,d} |".format(self.str_tag, "N=", total))
            ret.append(footerbar)
            ret.append("")
        else:
            ret = []
        if not self.generate_only:
            for line in ret:
                print(line)
        ret = "\n".join(ret)
        return ret


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
        self.bedoutzero = args.bedoutzero
        self.qlo = None
        self.qhi = None
        if args.outbeds != "outtab":
            self.qhi = args.qhi
            if args.qlo:
                try:
                    f = float(args.qlo)
                    self.qlo = f
                except Exception:
                    print("qlo not provided")
        nbw = len(args.bigwig)
        nlab = len(args.bigwiglabels)
        if nlab < nbw:
            self.bwlabels += ["Nolabel"] * (nbw - nlab)
        self.makeBed()

    def processVals(self, bw, isTop, isZero):
        """
        idea from http://gregoryzynda.com/python/numpy/contiguous/interval/2019/11/29/contiguous-regions.html
        Fast segmentation into regions by taking np.diff on the boolean array of over (under) cutpoint indicators in bwex.
        This only gives non-zero values at the segment boundaries where there's a change, so those zeros are all removed in bwexdnz
        leaving an array of segment start/end positions. That's twisted around into an array of start/end coordinates.
        Magical. Fast. Could do the same for means or medians over windows for sparse bigwigs like repeat regions.
        """
        if isTop:
            bwex = np.r_[False, bw >= self.bwtop, False]  # extend with 0s
        elif isZero:
            bwex = np.r_[False, bw == 0, False]  # extend with 0s
        else:
            bwex = np.r_[False, bw <= self.bwbot, False]
        bwexd = np.diff(bwex)
        bwexdnz = bwexd.nonzero()[0]  # start and end transition of each segment - nice!
        bwregions = np.reshape(bwexdnz, (-1, 2))
        return bwregions

    def writeBed(self, bed, bedfname):
        """
        potentially multiple
        """
        bed.sort()
        with open(bedfname, "w") as bedf:
            for b in bed:
                bedf.write("%s\t%d\t%d\t%s\t%d\n" % b)

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
        bedzero = []
        restab = []
        bwlabels = self.bwlabels
        bwnames = self.bwnames
        reshead = "bigwig\tcontig\tn\tmean\tstd\tmin\tmax\tqtop\tqbot"
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
            for chr in chrs:
                first_few = None
                bw = bwf.values(chr)
                values, counts = np.unique(bw, return_counts=True)
                nvalues = len(values)
                if nvalues <= 20:
                    histo = "\n".join(
                        [
                            "%s: %f occurs %d times" % (chr, values[x], counts[x])
                            for x in range(len(values))
                        ]
                    )
                else:
                    last10 = range(nvalues - 10, nvalues)
                    first_few = ["%.2f\t%d" % (values[x], counts[x]) for x in range(10)]
                    first_few += ["%.2f\t%d" % (values[x], counts[x]) for x in last10]
                    first_few.insert(0, "First/Last 10 value counts\nValue\tCount")
                    ha = asciihist(data=bw, bins=20, str_tag="%s_%s" % (bwlabel, chr))
                    histo = ha.draw()
                    histo = (
                        "\n".join(first_few)
                        + "\nHistogram of %s bigwig values\n" % bwlabel
                        + histo
                    )
                bw = bw[~np.isnan(bw)]  # some have NaN if parts of a contig not covered
                if self.bedoutzero is not None:
                    bwzero = self.processVals(bw, isTop=False, isZero=True)
                    for j, seg in enumerate(bwzero):
                        seglen = seg[1] - seg[0]
                        if seglen >= self.bedwin:
                            score = seglen
                            bedzero.append(
                                (
                                    chr,
                                    seg[0],
                                    seg[1],
                                    "%s_%d" % (bwlabel, score),
                                    score,
                                )
                            )
                if self.qhi is not None:
                    self.bwtop = np.quantile(bw, self.qhi)
                    bwhi = self.processVals(bw, isTop=True, isZero=False)
                    for j, seg in enumerate(bwhi):
                        seglen = seg[1] - seg[0]
                        if seglen >= self.bedwin:
                            score = np.sum(bw[seg[0]:seg[1]]) / float(seglen)
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
                    bwlo = self.processVals(bw, isTop=False, isZero=False)
                    for j, seg in enumerate(bwlo):
                        seglen = seg[1] - seg[0]
                        if seg[1] - seg[0] >= self.bedwin:
                            score = (
                                -1 * np.sum(bw[seg[0]:seg[1]]) / float(seglen)
                            )
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
                    resheadl = reshead.split("\t")
                    rowl = row.split()
                    desc = ["%s\t%s" % (resheadl[x], rowl[x]) for x in range(len(rowl))]
                    desc.insert(0, "Descriptive measures")
                    descn = "\n".join(desc)
                    restab.append(descn)
                    restab.append(histo)
        if os.path.isfile(fakepath):
            os.remove(fakepath)
        if self.tableoutfile:
            stable = "\n".join(restab)
            with open(self.tableoutfile, "w") as t:
                t.write(stable)
                t.write("\n")
        some = False
        if self.outbeds in ["outzero"]:
            self.writeBed(bedzero, self.bedoutzero)
            some = True
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
        if not ((self.outbeds == "outtab") or some):
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
    a("--bedoutzero", default=None)
    a("-w", "--bigwig", nargs="+")
    a("-n", "--bigwiglabels", nargs="+")
    a("-o", "--outbeds", default="outhilo", help="optional high and low combined bed")
    a("-t", "--tableoutfile", default=None)
    args = parser.parse_args()
    findOut(args)
