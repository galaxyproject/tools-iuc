import os

import ngsutils.support.ngs_utils
import pysam


class BedStreamer(object):
    '''
    Streams BedRegions from a BED file

    Note - this can only be used once! There is no mechanism to seek the stream.
    '''

    def __init__(self, fname=None, fileobj=None, quiet=False):
        if not fname and not fileobj:
            raise ValueError("You must specify either fname or fileobj!")

        self.reader = ngsutils.support.gzip_reader(fname=fname, quiet=quiet, fileobj=fileobj)

    def __iter__(self):
        return self

    def next(self):
        try:
            while True:
                line = self.reader.next().strip()
                if line and line[0] != '#':
                    cols = line.split('\t')
                    while len(cols) < 6:
                        cols.append('')

                    return BedRegion(*cols)
        except:
            raise StopIteration


class BedFile(object):
    '''
    BED files are read in their entirety memory, in a series of bins. Each bin
    is ~100kb in size. Each bin can then be iterated over.

    This is less efficient than using a proper index, but in reality, this
    usually isn't an issue. However, if the BED file has been Tabix indexed,
    that index will be used for random access.

    NOTE: This isn't very efficient, so perhaps this can be remodeled into a BedFile
    and a BedFileIndex where the file is indexed only if random access is requested.
    '''

    _bin_const = 100000

    def __init__(self, fname=None, fileobj=None, region=None):
        self._bins = {}
        self._bin_list = []
        self._cur_bin_idx = 0
        self._cur_bin_pos = 0
        self._tellpos = 0
        self._total = 0
        self._length = 0
        self.__tabix = None

        self.filename = fname

        if os.path.exists('%s.tbi' % fname):
            self.__tabix = pysam.Tabixfile(fname)

        if fileobj:
            self.__readfile(fileobj)
        elif fname:
            with ngsutils.support.ngs_utils.gzip_opener(fname) as fobj:
                self.__readfile(fobj)
        elif region:
            chrom, startend = region.split(':')
            if '-' in startend:
                start, end = [int(x) for x in startend.split('-')]
            else:
                start = int(startend)
                end = start
            start -= 1

            self.__add_region(BedRegion(chrom, start, end))
        else:
            raise ValueError("Must specify either filename, fileobj, or region")

    def __readfile(self, fobj):
        for line in fobj:
            line = line.strip()
            if line and line[0] != '#':
                cols = line.split('\t')
                while len(cols) < 6:
                    cols.append('')

                region = BedRegion(*cols)
                self.__add_region(region)

        self._bin_list.sort()
        for bin in self._bins:
            self._bins[bin].sort()

    def __add_region(self, region):
        self._total += region.end - region.start
        self._length += 1

        startbin = region.start / BedFile._bin_const
        endbin = region.end / BedFile._bin_const

        for bin in xrange(startbin, endbin + 1):
            if not (region.chrom, bin) in self._bins:
                self._bin_list.append((region.chrom, bin))
                self._bins[(region.chrom, bin)] = []
            self._bins[(region.chrom, bin)].append(region)

    def fetch(self, chrom, start, end, strand=None):
        '''
        For TABIX indexed BED files, find all regions w/in a range

        For non-TABIX index BED files, use the calculated bins, and
        output matching regions
        '''

        if self.__tabix:
            for match in self.__tabix.fetch(chrom, start, end):
                region = BedRegion(*match.split('\t'))
                if not strand or (strand and region.strand == strand):
                    yield region
        else:
            startbin = start / BedFile._bin_const
            endbin = end / BedFile._bin_const

            buf = set()

            for bin in xrange(startbin, endbin + 1):
                if (chrom, bin) in self._bins:
                    for region in self._bins[(chrom, bin)]:
                        if strand and strand != region.strand:
                            continue
                        if start <= region.start <= end or start <= region.end <= end:
                            if region not in buf:
                                yield region
                                buf.add(region)
                        elif region.start < start and region.end > end:
                            if region not in buf:
                                yield region
                                buf.add(region)

    def tell(self):
        return self._tellpos

    def close(self):
        pass

    @property
    def length(self):
        return self._length

    @property
    def total(self):
        return self._total

    def __iter__(self):
        self._cur_bin_idx = 0
        self._cur_bin_pos = 0
        self._tellpos = 0
        return self

    def next(self):
        if self._cur_bin_idx >= len(self._bin_list):
            raise StopIteration

        binvals = self._bins[self._bin_list[self._cur_bin_idx]]
        while self._cur_bin_pos < len(binvals):
            val = binvals[self._cur_bin_pos]
            self._cur_bin_pos += 1

            startbin = (val.chrom, val.start / BedFile._bin_const)
            if startbin == self._bin_list[self._cur_bin_idx]:
                self._tellpos += 1
                return val

        self._cur_bin_idx += 1
        self._cur_bin_pos = 0
        return self.next()


class BedRegion(object):
    def __init__(self, chrom, start, end, name='', score='', strand='', thickStart='', thickEnd='', rgb='', *args):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name

        if score == '':
            self.score = 0
        else:
            self.score = float(score)

        if strand == '':
            self.strand = None
        else:
            self.strand = strand

        if thickStart == '':
            self.thickStart = None
        else:
            self.thickStart = thickStart

        if thickEnd == '':
            self.thickEnd = None
        else:
            self.thickEnd = thickEnd

        if rgb == '':
            self.rgb = None
        else:
            self.rgb = rgb

        self.extras = args

    def clone(self, chrom=None, start=None, end=None, name=None, score=None, strand=None, thickStart=None, thickEnd=None, rgb=None, *args):
        cols = []
        cols.append(self.chrom if chrom is None else chrom)
        cols.append(self.start if start is None else start)
        cols.append(self.end if end is None else end)
        cols.append(self.name if name is None else name)
        cols.append(self.score if score is None else score)
        cols.append(self.strand if strand is None else strand)
        cols.append(self.thickStart if thickStart is None else thickStart)
        cols.append(self.thickEnd if thickEnd is None else thickEnd)
        cols.append(self.rgb if rgb is None else rgb)

        for i, val in enumerate(self.extras):
            if len(args) > i:
                cols.append(args[i])
            else:
                cols.append(val)

        return BedRegion(*cols)

    @property
    def score_int(self):
        score = str(self.score)
        if score[-2:] == '.0':
            score = score[:-2]

        return score

    def __key(self):
        return (self.chrom, self.start, self.end, self.strand, self.name)

    def __lt__(self, other):
        return self.__key() < other.__key()

    def __gt__(self, other):
        return self.__key() > other.__key()

    def __eq__(self, other):
        return self.__key() == other.__key()

    def write(self, out):
        out.write('%s\n' % self)

    def __repr__(self):
        outcols = []

        if self.rgb:
            outcols.append(self.rgb)
        if self.thickEnd or outcols:
            outcols.append(self.thickEnd if self.thickEnd else self.end)
        if self.thickStart or outcols:
            outcols.append(self.thickStart if self.thickStart else self.start)
        if self.strand or outcols:
            outcols.append(self.strand)
        if self.score_int != '' or outcols:
            outcols.append(self.score_int)
        if self.name or outcols:
            outcols.append(self.name)

        outcols.append(self.end)
        outcols.append(self.start)
        outcols.append(self.chrom)

        return '\t'. join([str(x) for x in outcols[::-1]])
