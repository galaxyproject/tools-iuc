class RangeMatch(object):
    '''
    Simple genomic ranges.  You can define chrom:start-end ranges, then ask if a
    particular genomic coordinate maps to any of those ranges.  This is less-
    efficient than an R-Tree, but easier to code.
    '''

    def __init__(self, name):
        self.ranges = {}
        self.name = name

    def add_range(self, chrom, strand, start, end):
        if chrom not in self.ranges:
            self.ranges[chrom] = {}

        bin = start / 100000
        if bin not in self.ranges[chrom]:
            self.ranges[chrom][bin] = []
        self.ranges[chrom][bin].insert(0, (start, end, strand))

        if (end / 100000) != bin:
            for bin in xrange(bin + 1, (end / 100000) + 1):
                if bin not in self.ranges[chrom]:
                    self.ranges[chrom][bin] = []
                self.ranges[chrom][bin].insert(0, (start, end, strand))

    def get_tag(self, chrom, strand, pos, ignore_strand=False):
        '''
        returns (region, is_reverse_orientation)
        '''
        if chrom not in self.ranges:
            return None, False
        bin = pos / 100000
        if bin not in self.ranges[chrom]:
            return None, False
        for start, end, r_strand in self.ranges[chrom][bin]:
            if pos >= start and pos <= end:
                if ignore_strand or strand == r_strand:
                    return self.name, False
                return self.name, True
        return None, False


class RegionTagger(object):
    def __init__(self, gtf, valid_chroms=None, only_first_fragment=True):
        self.regions = []
        self.counts = {}
        self.only_first_fragment = only_first_fragment

        coding = RangeMatch('coding')
        exons = RangeMatch('other-exon')
        utr_5 = RangeMatch('utr-5')
        utr_3 = RangeMatch('utr-3')
        introns = RangeMatch('intron')
        promoters = RangeMatch('promoter')

        for gene in gtf.genes:
            if valid_chroms and gene.chrom not in valid_chroms:
                continue
            if gene.strand == '+':
                promoters.add_range(gene.chrom, gene.strand, gene.start - 2000, gene.start)
            else:
                promoters.add_range(gene.chrom, gene.strand, gene.end, gene.end + 2000)

            for transcript in gene.transcripts:
                if transcript.has_cds:
                    for start, end in transcript.cds:
                        coding.add_range(gene.chrom, gene.strand, start, end)

                    # TODO: Fix this so that it iterates over exons in the 5'/3' UTRS
                    for s, e in transcript.utr_5:
                        utr_5.add_range(gene.chrom, gene.strand, s, e)
                    for s, e in transcript.utr_3:
                        utr_3.add_range(gene.chrom, gene.strand, s, e)

                last_end = None
                for start, end in transcript.exons:
                    if last_end:
                        introns.add_range(gene.chrom, gene.strand, last_end, start)
                    exons.add_range(gene.chrom, gene.strand, start, end)
                    last_end = end

        self.regions.append(coding)
        self.regions.append(utr_5)
        self.regions.append(utr_3)
        self.regions.append(exons)
        self.regions.append(introns)
        self.regions.append(promoters)

        self.counts['coding'] = 0
        self.counts['coding-rev'] = 0
        self.counts['other-exon'] = 0
        self.counts['utr-5'] = 0
        self.counts['utr-3'] = 0
        self.counts['utr-5-rev'] = 0
        self.counts['utr-3-rev'] = 0
        self.counts['intron'] = 0
        self.counts['promoter'] = 0
        self.counts['other-exon-rev'] = 0
        self.counts['intron-rev'] = 0
        self.counts['promoter-rev'] = 0
        self.counts['junction'] = 0
        self.counts['intergenic'] = 0
        self.counts['mitochondrial'] = 0

    def add_read(self, read, chrom):
        if read.is_unmapped:
            return

        if self.only_first_fragment and read.is_paired and not read.is_read1:
            return

        tag = None
        is_rev = False

        strand = '-' if read.is_reverse else '+'

        if chrom == 'chrM':
            tag = 'mitochondrial'

        if not tag:
            for op, length in read.cigar:
                if op == 3:
                    tag = 'junction'
                    break

        if not tag:
            for region in self.regions:
                tag, is_rev = region.get_tag(chrom, strand, read.pos)
                if tag:
                    break

        if not tag:
            tag = 'intergenic'

        if tag:
            if is_rev:
                self.counts['%s-rev' % tag] += 1
            else:
                self.counts[tag] += 1

        return tag

    def tag_region(self, chrom, start, end, strand):
        tag = None
        is_rev = False

        if chrom == 'chrM' or chrom == 'M':
            tag = 'mitochondrial'

        if not tag:
            for region in self.regions:
                tag, is_rev = region.get_tag(chrom, strand, start)
                if is_rev:
                    tag = '%s-rev' % tag

                if start != end:
                    endtag, is_rev = region.get_tag(chrom, strand, end)
                    if is_rev:
                        endtag = '%s-rev' % endtag

                if tag and endtag and endtag != tag:
                    tag = '%s/%s' % (tag, endtag)

        if not tag:
            tag = 'intergenic'
