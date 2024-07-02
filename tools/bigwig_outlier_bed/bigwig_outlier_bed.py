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
import numpy as np
import pybigtools
import sys


from pathlib import Path

class findOut():

    def __init__(self, args):
        self.bwnames=args.bigwig
        self.bwlabels=args.bigwiglabels
        self.bedwin=args.minwin
        self.qlo=args.qlo
        self.qhi=args.qhi
        self.bedouthilo=args.bedouthilo
        self.bedouthi=args.bedouthi
        self.bedoutlo=args.bedoutlo
        self.tableout = args.tableout
        self.bedwin = args.minwin
        self.qhi = args.qhi
        self.qlo = args.qlo
        self.makeBed()

    def processVals(self, bw, isTop):
        # http://gregoryzynda.com/python/numpy/contiguous/interval/2019/11/29/contiguous-regions.html
        if isTop:
            bwex = np.r_[False, bw >= self.bwtop, False] # extend with 0s
        else:
            bwex = np.r_[False, bw <= self.bwbot, False]
        bwexd = np.diff(bwex)
        bwexdnz = bwexd.nonzero()[0]
        bwregions = np.reshape(bwexdnz, (-1,2))
        return bwregions

    def writeBed(self, bed, bedfname):
        """
        potentially multiple
        """
        bed.sort()
        beds = ['%s\t%d\t%d\t%s\t%d' % x for x in bed]
        with open(bedfname, "w") as bedf:
            bedf.write('\n'.join(beds))
            bedf.write('\n')
        print('Wrote %d bed regions to %s' % (len(bed), bedfname))
        
    def makeBed(self):
        bedhi = []
        bedlo = []
        bwlabels = self.bwlabels
        bwnames = self.bwnames
        print('bwnames=', bwnames, "bwlabs=", bwlabels)
        for i, bwname in enumerate(bwnames):
            bwlabel = bwlabels[i].replace(" ",'')
            p = Path('in.bw')
            p.symlink_to( bwname ) # required by pybigtools (!)
            bwf = pybigtools.open('in.bw')
            chrlist = bwf.chroms()
            chrs = list(chrlist.keys())
            chrs.sort()
            restab = ["contig\tn\tmean\tstd\tmin\tmax\tqtop\tqbot"]
            for chr in chrs:
                bw = bwf.values(chr)
                bw = bw[~np.isnan(bw)] # some have NaN if parts of a contig not covered
                if self.qhi is not None:
                    self.bwtop = np.quantile(bw, self.qhi)
                    bwhi = self.processVals(bw, isTop=True)
                    for i, seg in enumerate(bwhi):
                        if seg[1] - seg[0] >= self.bedwin:
                            bedhi.append((chr, seg[0], seg[1], '%s_hi' % (bwlabel), 1))
                if self.qlo is not None:
                    self.bwbot = np.quantile(bw, self.qlo)
                    bwlo = self.processVals(bw, isTop=False)            
                    for i, seg in enumerate(bwlo):
                        if seg[1] - seg[0] >= self.bedwin:
                            bedlo.append((chr, seg[0], seg[1], '%s_lo' % (bwlabel), -1))
                bwmean = np.mean(bw)
                bwstd = np.std(bw)
                bwmax = np.max(bw)
                nrow = np.size(bw)
                bwmin = np.min(bw)
                restab.append('%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f' % (chr,nrow,bwmean,bwstd,bwmin,bwmax,self.bwtop,self.bwbot))        
        print('\n'.join(restab), '\n')
        if self.tableout:
            with open(self.tableout) as t:
                t.write('\n'.join(restab))
                t.write('\n')
        if self.bedoutlo:
            if self.qlo:
                self.writeBed(bedlo, self.bedoutlo)
        if self.bedouthi:
            if self.qhi:
                self.writeBed(bedhi, self.bedouthi)
        if self.bedouthilo:
            allbed = bedlo + bedhi
            self.writeBed(allbed, self.bedouthilo)
        return restab
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    a = parser.add_argument
    a('-m', '--minwin',default=10, type=int)
    a('-l', '--qlo',default=None, type=float)
    a('-i', '--qhi',default=None, type=float)
    a('-w', '--bigwig', nargs='+')
    a('-n', '--bigwiglabels', nargs='+')
    a('-o', '--bedouthilo', default=None, help="optional high and low combined bed")
    a('-u', '--bedouthi', default=None, help="optional high only bed")
    a('-b', '--bedoutlo', default=None, help="optional low only bed")
    a('-t', '--tableout', default=None)
    args = parser.parse_args()
    print('args=', args)
    if not (args.bedouthilo or args.bedouthi or args.bedoutlo):
        sys.stderr.write("bigwig_outlier_bed.py cannot usefully run - need a bed output choice - must be one of low only, high only or both combined")
        sys.exit(2)
    if not (args.qlo or args.qhi):
        sys.stderr.write("bigwig_outlier_bed.py cannot usefully run - need one or both of quantile cutpoints qhi and qlo")
        sys.exit(2)
    restab = findOut(args)
    if args.tableout:
        with open(args.tableout, 'w') as tout:
            tout.write('\n'.join(restab))
            tout.write('\n')
    

