#!/usr/bin/env python
# category General
# desc Removes reads from a BAM file based on criteria
"""
Removes reads from a BAM file based on criteria

Given a BAM file, this script will only allow reads that meet filtering
criteria to be written to output. The output is another BAM file with the
reads not matching the criteria removed.

Note: this does not adjust tag values reflecting any filtering. (for example:
      if a read mapped to two locations (IH:i:2), and one was removed by
      filtering, the IH:i tag would still read IH:i:2).

Currently, the available filters are:
    -minlen val                Remove reads that are smaller than {val}
    -maxlen val                Remove reads that are larger than {val}
    -mapped                    Keep only mapped reads
    -unmapped                  Keep only unmapped reads
    -properpair                Keep only properly paired reads (both mapped,
                               correct orientation, flag set in BAM)
    -noproperpair              Keep only not-properly paired reads

    -mask bitmask              Remove reads that match the mask (base 10/hex)
    -uniq {length}             Remove reads that are have the same sequence
                               Note: BAM file should be sorted
                               (up to an optional length)
    -uniq_start                Remove reads that start at the same position
                               Note: BAM file should be sorted
                               (Use only for low-coverage samples)

    -mismatch num              # mismatches or indels
                               indel always counts as 1 regardless of length
                               (requires NM tag in reads)

    -mismatch_dbsnp num dbsnp.txt.bgz
                               # mismatches or indels - not in dbSNP.
                               Variations are called using the MD tag.
                               Variations that are found in the dbSNP list are
                               not counted as mismatches. The dbSNP list is a
                               Tabix-indexed dump of dbSNP (from UCSC Genome
                               Browser). Indels in dbSNP are also counted.
                               Adds a 'ZS:i' tag with the number of found SNPs
                               in the read.
                               (requires NM and MD tags)

                               Example command for indexing:
                               ngsutils tabixindex snp.txt.gz -s 2 -b 3 -e 4 -0

    -mismatch_ref num ref.fa   # mismatches or indel - looks up mismatches
                               directly in a reference FASTA file
                               (use if NM tag not present)

    -mismatch_ref_dbsnp num ref.fa dbsnp.txt.bgz
                               # mismatches or indels - looks up mismatches
                               directly from a reference FASTA file. (See
                               -mismatch_dbsnp for dbSNP matching)
                               (use if NM or MD tag not present)

    -nosecondary               Remove reads that have the 0x100 flag set
    -noqcfail                  Remove reads that have the 0x200 flag set
    -nopcrdup                  Remove reads that have the 0x400 flag set


    -exclude ref:start-end     Remove reads in this region (1-based start)
    -excludebed file.bed {nostrand}
                               Remove reads that are in any of the regions
                               from the given BED file. If 'nostrand' is given,
                               strand information from the BED file is ignored.

    -include ref:start-end     Remove reads NOT in the region (can only be one)
    -includebed file.bed {nostrand}
                               Remove reads that are NOT any of the regions
                               from the given BED file. If 'nostrand' is given,
                               strand information from the BED file is ignored.

                               Note: If this is a large dataset, use
                               "bamutils extract" instead.

    -includeref refname        Exclude reads NOT mapped to a reference
    -excluderef refname        Exclude reads mapped to a particular reference
                               (e.g. chrM, or _dup chromosomes)

    -whitelist fname           Remove reads that aren't on this list (by name)
    -blacklist fname           Remove reads that are on this list (by name)
                                 These lists can be whitespace-delimited with
                                 the read name as the first column.
    -maximum_mismatch_ratio val
                               Filter by maximum mismatch ratio (fraction of length)

    -eq  tag_name value
    -lt  tag_name value
    -lte tag_name value
    -gt  tag_name value
    -gte tag_name value

    As a special case, "MAPQ" can be used as the tag_name and the SAM MAPQ
    value will be used.

Common tags to filter by:
    AS      Alignment score
    IH      Number of alignments
    NM      Edit distance (each indel counts as many as its length)

    MAPQ    Mapping quality (defined as part of SAM spec)

    The tag type (:i, :f, :Z) is optional.

"""
from __future__ import print_function

import os
import sys

import pysam
from ngsutils.bam import bam_iter, read_calc_mismatches, read_calc_mismatches_gen, read_calc_mismatches_ref, read_calc_variations
from ngsutils.bed import BedFile
from ngsutils.support.dbsnp import DBSNP


def usage():
    print(__doc__)
    print("""
Usage: bamutils filter in.bam out.bam {-failed out.txt} criteria...

Options:
  -failed fname    A text file containing the read names of all reads
                   that were removed with filtering

Example:
bamutils filter filename.bam output.bam -mapped -gte AS:i 1000

This will remove all unmapped reads, as well as any reads that have an AS:i
value less than 1000.
""")
    sys.exit(1)


class Unique(object):
    def __init__(self, length=None):
        if length:
            self.length = int(length)
        else:
            self.length = None

        self.last_pos = None
        self.pos_reads = set()

    def __repr__(self):
        return "uniq"

    def filter(self, bam, read):
        if self.last_pos != (read.tid, read.pos):
            self.last_pos = (read.tid, read.pos)
            self.pos_reads = set()

        if read.is_reverse:
            seq = read.seq[::-1]  # ignore revcomp for now, it isn't needed, just need to compare in the proper order
        else:
            seq = read.seq

        if self.length:
            seq = seq[:self.length]

        if seq in self.pos_reads:
            return False

        self.pos_reads.add(seq)
        return True

    def close(self):
        pass


class UniqueStart(object):
    def __init__(self):
        self.last_tid = None
        self.last_fwd_pos = -1
        self.rev_pos_list = None

    def __repr__(self):
        return "uniq_start"

    def filter(self, bam, read):
        if self.last_tid is None or self.last_tid != read.tid:
            self.last_tid = read.tid
            self.rev_pos = set()
            self.last_fwd_pos = -1
            self.last_rev_pos = -1

        if read.is_reverse:
            # check reverse reads from their start (3' aend)
            # these aren't necessarily in the correct
            # order in the file, so we have to track them in a set

            start_pos = read.aend

            # clean up hash if necesary
            # Allow up to 100k over, to balance cleaning up the rev_pos hash
            # and memory
            if read.pos > (self.last_rev_pos + 100000):
                self.last_rev_pos = start_pos
                del_list = []
                for k in self.rev_pos:
                    if k < read.pos:
                        del_list.append(k)

                for k in del_list:
                    self.rev_pos.remove(k)

            if start_pos not in self.rev_pos:
                self.rev_pos.add(start_pos)
                return True
            return False
        else:
            if read.pos != self.last_fwd_pos:
                self.last_fwd_pos = read.pos
                return True
            return False

        # if self.last_pos != (read.tid, read.pos):
        #     self.last_pos = (read.tid, read.pos)
        #     return True
        # return False

    def close(self):
        pass


class Blacklist(object):
    def __init__(self, fname):
        self.fname = fname
        self.notallowed = set()
        with open(fname) as f:
            for line in f:
                self.notallowed.add(line.strip().split()[0])

    def filter(self, bam, read):
        return read.qname not in self.notallowed

    def __repr__(self):
        return 'Blacklist: %s' % (self.fname)

    def close(self):
        pass


class Whitelist(object):
    def __init__(self, fname):
        self.fname = fname
        self.allowed = set()
        with open(fname) as f:
            for line in f:
                self.allowed.add(line.strip().split()[0])

    def filter(self, bam, read):
        return read.qname in self.allowed

    def __repr__(self):
        return 'Whitelist: %s' % (self.fname)

    def close(self):
        pass


class IncludeRegion(object):
    _excludes = []
    _last = None

    def __init__(self, region):
        IncludeRegion._excludes.append(ExcludeRegion(region))
        # self.excl = ExcludeRegion(region)

    def filter(self, bam, read):
        if read == IncludeRegion._last:
            return True

        IncludeRegion._last = read

        for excl in IncludeRegion._excludes:
            if not excl.filter(bam, read):
                return True
        return False

        # return not self.excl.filter(bam,read)
    def __repr__(self):
        return 'Including: %s' % (self.excl.region)

    def close(self):
        pass


class IncludeBED(object):
    def __init__(self, fname, nostrand=None):
        self.excl = ExcludeBED(fname, nostrand)

    def filter(self, bam, read):
        return not self.excl.filter(bam, read)

    def __repr__(self):
        return 'Including from BED: %s%s' % (self.excl.fname, ' nostrand' if self.excl.nostrand else '')

    def close(self):
        pass


class ExcludeRegion(object):
    def __init__(self, region):
        self.region = region
        spl = region.split(':')
        self.chrom = spl[0]
        se = [int(x) for x in spl[1].split('-')]
        self.start = se[0] - 1
        self.end = se[1]

    def filter(self, bam, read):
        if not read.is_unmapped:
            if bam.getrname(read.tid) == self.chrom:
                if self.start <= read.pos <= self.end:
                    return False
                if self.start <= read.aend <= self.end:
                    return False
        return True

    def __repr__(self):
        return 'Excluding: %s' % (self.region)

    def close(self):
        pass


class ExcludeRef(object):
    def __init__(self, ref):
        self.ref = ref

    def filter(self, bam, read):
        if not read.is_unmapped:
            if bam.getrname(read.tid) == self.ref:
                    return False
        return True

    def __repr__(self):
        return 'Excluding: %s' % (self.ref)

    def close(self):
        pass


class IncludeRef(object):
    def __init__(self, ref):
        self.ref = ref

    def filter(self, bam, read):
        if not read.is_unmapped:
            if bam.getrname(read.tid) == self.ref:
                    return True
        return False

    def __repr__(self):
        return 'Including: %s' % (self.ref)

    def close(self):
        pass


class ExcludeBED(object):
    def __init__(self, fname, nostrand=None):
        self.regions = {}  # store BED regions as keyed bins (chrom, bin)
        self.fname = fname
        if nostrand == 'nostrand':
            self.nostrand = True
        else:
            self.nostrand = False

        self.bed = BedFile(fname)
        # with open(fname) as f:
        #     for line in f:
        #         if not line:
        #             continue
        #         if line[0] == '#':
        #             continue
        #         cols = line.strip().split('\t')

        #         chrom = cols[0]
        #         start = int(cols[1])
        #         end = int(cols[2])
        #         if self.nostrand:
        #             strand = '?'
        #         else:
        #             strand = cols[5]

        #         startbin = start / 100000
        #         endbin = end / 100000

        #         for bin in xrange(startbin, endbin + 1):
        #             if not (chrom, bin) in self.regions:
        #                 self.regions[(chrom, bin)] = []
        #         self.regions[(chrom, bin)].append((start, end, strand))

    def filter(self, bam, read):
        if not read.is_unmapped:
            if self.nostrand:
                strand = None
            elif read.is_reverse:
                strand = '-'
            else:
                strand = '+'

            for region in self.bed.fetch(bam.getrname(read.tid), read.pos, read.aend, strand):
                # region found, exclude read
                return False
            return True

        #     bin = read.pos / 100000
        #     ref = bam.getrname(read.tid)

        #     if not (ref, bin) in self.regions:
        #         return True

        #     for start, end, strand in self.regions[(ref, bin)]:
        #         if not self.nostrand:
        #             if strand == '+' and read.is_reverse:
        #                 continue
        #             if strand == '-' and not read.is_reverse:
        #                 continue
        #         if start <= read.pos <= end:
        #             return False
        #         if start <= read.aend <= end:
        #             return False
        # return True

    def __repr__(self):
        return 'Excluding from BED: %s%s' % (self.fname, ' nostrand' if self.nostrand else '')

    def close(self):
        pass


class Mismatch(object):
    def __init__(self, num):
        self.num = int(num)

    def filter(self, bam, read):
        if read.is_unmapped:
            return False

        if read_calc_mismatches(read) > self.num:
            return False

        return True

    def __repr__(self):
        return '>%s mismatch%s' % (self.num, '' if self.num == 1 else 'es')

    def close(self):
        pass


class MismatchRef(object):
    def __init__(self, num, refname):
        self.num = int(num)
        self.refname = refname

        if not os.path.exists('%s.fai' % refname):
            pysam.faidx(refname)

        self.ref = pysam.Fastafile(refname)

    def filter(self, bam, read):
        if read.is_unmapped:
            return False

        chrom = bam.getrname(read.tid)
        if read_calc_mismatches_ref(self.ref, read, chrom) > self.num:
            return False

        return True

    def __repr__(self):
        return '>%s mismatch%s in %s' % (self.num, '' if self.num == 1 else 'es', os.path.basename(self.refname))

    def close(self):
        self.ref.close()


class MismatchDbSNP(object):
    def __init__(self, num, fname, verbose=None):
        sys.stderr.write('Note: MismatchDbSNP is considered *experimental*\n')

        self.num = int(num)
        self.fname = fname
        self.dbsnp = DBSNP(fname)
        if verbose == 'verbose':
            self.verbose = True
        else:
            self.verbose = False

    def filter(self, bam, read):
        if read.is_unmapped:
            return False

        if read_calc_mismatches(read) <= self.num:
            return True

        chrom = bam.getrname(read.tid)

        mm = 0
        snps = 0

        for op, pos, seq in read_calc_variations(read):
            if not self.dbsnp.is_valid_variation(chrom, op, pos, seq, self.verbose):
                mm += 1
            else:
                snps += 1

        if mm > self.num:
            return False

        if snps:
            read.tags = read.tags + [('ZS', snps)]

        return True

    def __repr__(self):
        return '>%s mismatch%s using %s' % (self.num, '' if self.num == 1 else 'es', os.path.basename(self.fname))

    def close(self):
        self.dbsnp.close()


class MismatchRefDbSNP(object):
    def __init__(self, num, refname, dbsnpname):
        sys.stderr.write('Note: MismatchRefDbSNP is considered *experimental*\n')
        self.num = int(num)
        self.refname = refname
        self.dbsnp = DBSNP(dbsnpname)

        if not os.path.exists('%s.fai' % refname):
            pysam.faidx(refname)

        self.ref = pysam.Fastafile(refname)

    def filter(self, bam, read):
        if read.is_unmapped:
            return False

        chrom = bam.getrname(read.tid)

        mm = 0
        snps = 0

        for op, pos, seq in read_calc_mismatches_gen(self.ref, read, chrom):
            if not self.dbsnp.is_valid_variation(chrom, op, pos, seq):
                mm += 1
            else:
                snps += 1

        if mm > self.num:
            return False

        if snps:
            read.tags = read.tags + [('ZS', snps)]

        return True

    def __repr__(self):
        return '>%s mismatch%s using %s/%s' % (self.num, '' if self.num == 1 else 'es', os.path.basename(self.dbsnpname), os.path.basename(self.refname))

    def close(self):
        self.ref.close()
        self.dbsnp.close()


class Mapped(object):
    def __init__(self):
        pass

    def filter(self, bam, read):
        if read.is_paired and (read.is_unmapped or read.mate_is_unmapped):
            return False
        elif read.is_unmapped:
            return False
        return True

    def __repr__(self):
        return 'is mapped'

    def close(self):
        pass


class Unmapped(object):
    def __init__(self):
        pass

    def filter(self, bam, read):
        if read.is_paired and (read.is_unmapped or read.mate_is_unmapped):
            return True
        elif read.is_unmapped:
            return True
        return False

    def __repr__(self):
        return 'is unmapped'

    def close(self):
        pass


class ProperPair(object):
    def __init__(self):
        pass

    def filter(self, bam, read):
        if not read.is_paired:
            return False

        if read.is_unmapped or read.mate_is_unmapped:
            return False

        if read.is_reverse == read.mate_is_reverse:
            return False

        return read.is_proper_pair

    def __repr__(self):
        return 'proper pair'

    def close(self):
        pass


class NoProperPair(object):
    def __init__(self):
        self.proper = ProperPair()
        pass

    def filter(self, bam, read):
        return not self.proper.filter(bam, read)

    def __repr__(self):
        return 'not proper pairs'

    def close(self):
        pass


class MaskFlag(object):
    def __init__(self, value):
        if isinstance(value, int):
            self.flag = value
        else:
            if value[0:2] == '0x':
                self.flag = int(value, 16)
            else:
                self.flag = int(value)

    def __repr__(self):
        return "Doesn't match flag: %s" % self.flag

    def filter(self, bam, read):
        return (read.flag & self.flag) == 0

    def close(self):
        pass


class SecondaryFlag(object):
    def __repr__(self):
        return "no 0x100 (secondary) flag"

    def filter(self, bam, read):
        return not read.is_secondary

    def close(self):
        pass


class ReadMinLength(object):
    def __init__(self, minval):
        self.minval = int(minval)

    def __repr__(self):
        return "read length min: %s" % self.minval

    def filter(self, bam, read):
        return len(read.seq) >= self.minval

    def close(self):
        pass


class ReadMaxLength(object):
    def __init__(self, val):
        self.val = int(val)

    def __repr__(self):
        return "read length max: %s" % self.val

    def filter(self, bam, read):
        return len(read.seq) <= self.val

    def close(self):
        pass


class MaximumMismatchRatio(object):
    def __init__(self, ratio):
        self.ratio = float(ratio)

    def __repr__(self):
        return "maximum mismatch ratio: %s" % self.val

    def filter(self, bam, read):
        return read_calc_mismatches(read) <= self.ratio * len(read.seq)

    def close(self):
        pass


class QCFailFlag(object):
    def __repr__(self):
        return "no 0x200 (qcfail) flag"

    def filter(self, bam, read):
        return not read.is_qcfail

    def close(self):
        pass


class PCRDupFlag(object):
    def __repr__(self):
        return "no 0x400 (pcrdup) flag"

    def filter(self, bam, read):
        return not read.is_duplicate

    def close(self):
        pass


class _TagCompare(object):
    def __init__(self, tag, value):
        self.args = '%s %s' % (tag, value)

        if ':' in tag:
            self.tag = tag.split(':')[0]
            tagtype = tag.split(':')[1]
            if tagtype == 'i':
                self.value = int(value)
            elif tagtype == 'f':
                self.value = float(value)
            elif tagtype == 'H':
                self.value = int(value)  # not sure about this
            else:
                self.value = value
        else:
            self.tag = tag

            # guess at type...
            try:
                self.value = int(value)
            except ValueError:
                try:
                    self.value = float(value)
                except ValueError:
                    self.value = value

    def get_value(self, read):
        if self.tag == 'MAPQ':
            return read.mapq

        for name, value in read.tags:
            if name == self.tag:
                return value

        return None

    def __repr__(self):
        return "%s %s %s" % (self.tag, self.__class__.op, self.value)

    def close(self):
        pass


class TagLessThan(_TagCompare):
    op = '<'

    def filter(self, bam, read):
        if self.get_value(read) < self.value:
            return True
        return False


class TagLessThanEqual(_TagCompare):
    op = '<='

    def filter(self, bam, read):
        if self.get_value(read) <= self.value:
            return True
        return False


class TagGreaterThan(_TagCompare):
    op = '>'

    def filter(self, bam, read):
        if self.get_value(read) > self.value:
            return True
        return False


class TagGreaterThanEqual(_TagCompare):
    op = '>='

    def filter(self, bam, read):
        if self.get_value(read) >= self.value:
            return True
        return False


class TagEqual(_TagCompare):
    op = '='

    def filter(self, bam, read):
        if self.get_value(read) == self.value:
            return True
        return False


_criteria = {
    'mapped': Mapped,
    'unmapped': Unmapped,
    'properpair': ProperPair,
    'noproperpair': NoProperPair,
    'noqcfail': QCFailFlag,
    'nosecondary': SecondaryFlag,
    'nopcrdup': PCRDupFlag,
    'mask': MaskFlag,
    'lt': TagLessThan,
    'gt': TagGreaterThan,
    'lte': TagLessThanEqual,
    'gte': TagGreaterThanEqual,
    'eq': TagEqual,
    'mismatch': Mismatch,
    'mismatch_ref': MismatchRef,
    'mismatch_dbsnp': MismatchDbSNP,
    'mismatch_ref_dbsnp': MismatchRefDbSNP,
    'exclude': ExcludeRegion,
    'excludebed': ExcludeBED,
    'excluderef': ExcludeRef,
    'includeref': IncludeRef,
    'include': IncludeRegion,
    'includebed': IncludeBED,
    'whitelist': Whitelist,
    'blacklist': Blacklist,
    'length': ReadMinLength,  # deprecated
    'minlen': ReadMinLength,
    'maxlen': ReadMaxLength,
    'uniq': Unique,
    'maximum_mismatch_ratio': MaximumMismatchRatio
}


def bam_filter(infile, outfile, criteria, failedfile=None, verbose=False):
    if verbose:
        sys.stderr.write('Input file  : %s\n' % infile)
        sys.stderr.write('Output file : %s\n' % outfile)
        if failedfile:
            sys.stderr.write('Failed reads: %s\n' % failedfile)
        sys.stderr.write('Criteria:\n')
        for criterion in criteria:
            sys.stderr.write('    %s\n' % criterion)

        sys.stderr.write('\n')

    bamfile = pysam.Samfile(infile, "rb")
    outfile = pysam.Samfile(outfile, "wb", template=bamfile)

    if failedfile:
        failed_out = open(failedfile, 'w')
    else:
        failed_out = None

    passed = 0
    failed = 0

    def _callback(read):
        return "%s | %s kept,%s failed" % ('%s:%s' % (bamfile.getrname(read.tid), read.pos) if read.tid > -1 else 'unk', passed, failed)

    for read in bam_iter(bamfile, quiet=True):
        p = True

        for criterion in criteria:
            if not criterion.filter(bamfile, read):
                p = False
                failed += 1
                if failed_out:
                    failed_out.write('%s\t%s\n' % (read.qname, criterion))
                # outfile.write(read_to_unmapped(read))
                break
        if p:
            passed += 1
            outfile.write(read)

    bamfile.close()
    outfile.close()
    if failed_out:
        failed_out.close()
    sys.stdout.write("%s kept\n%s failed\n" % (passed, failed))

    for criterion in criteria:
        criterion.close()


def read_to_unmapped(read):
    r'''
    Example unmapped read

    2358_2045_1839    |  4 | *                             |  0    |  0   | *                  | *  |  0 |  0 | GGACGAGAAGGAGTATTTCTCCGAGAACACATTCACGGAGAGTCTAACTC           | 0QT\VNO^IIRJKXTIHIKTY\STZZ[XOJKPWLHJQQQ^XPQYIIRQII           | PG:Z:bfast   | AS:i:-                     | NH:i:0   | IH:i:1   | HI:i:1                                     | CS:Z:T10213222020221330022203222011113021130222212230122             | CQ:Z:019=2+><%@.'>52%)'15?90<7;:5+(-49('-7-<>5-@5%<.2%=            | CM:i:0   | XA:i:4

    Example mapped read:

    2216_1487_198     |  16 | chr11:19915291-19925392      |  5531 |   12 | 24M2D26M           | *  |  0 |  0 | TCTGTCTGGGTGCAGTAGCTATACGTAATCCCAGCACTTTGGGAGGCCAA           | 1C8FZ`M""""""WZ""%"#\I;"`R@D]``ZZ``\QKX\Y]````IK^`           | PG:Z:bfast   | AS:i:1150   | NM:i:2   | NH:i:10   | IH:i:1   | HI:i:1   | MD:Z:24^CT26                   | CS:Z:T00103022001002113210023031213331122121223321221122             | CQ:Z:A8=%;AB<;97<3/9:>>3>@?5&18@-%;866:94=14:=?,%?:7&)1            | CM:i:9   | XA:i:4   | XE:Z:-------23322---2-11----2----------------------------
    '''
    # TODO: rewrite read properties to unmapped state
    #       if colorspace: convert CS to NT directly
    #       remove excess tags

    read.is_unmapped = True
    read.tid = None
    read.pos = 0
    read.mapq = 0
    return read


if __name__ == '__main__':
    infile = None
    outfile = None
    failed = None
    criteria = []

    crit_args = []
    last = None
    verbose = False
    fail = False

    for arg in sys.argv[1:]:
        if last == '-failed':
            failed = arg
            last = None
        elif arg == '-h':
            usage()
        elif arg == '-failed':
            last = arg
        elif arg == '-v':
            verbose = True
        elif not infile and os.path.exists(arg):
            infile = arg
        elif not outfile:
            outfile = arg
        elif arg[0] == '-':
            if not arg[1:] in _criteria:
                print("Unknown criterion: %s" % arg)
                fail = True
            if crit_args:
                criteria.append(_criteria[crit_args[0][1:]](*crit_args[1:]))
            crit_args = [arg, ]
        elif crit_args:
            crit_args.append(arg)
        else:
            print("Unknown argument: %s" % arg)
            fail = True

    if not fail and crit_args:
        criteria.append(_criteria[crit_args[0][1:]](*crit_args[1:]))

    if fail or not infile or not outfile or not criteria:
        if not infile and not outfile and not criteria:
            usage()

        if not infile:
            print("Missing: input bamfile")
        if not outfile:
            print("Missing: output bamfile")
        if not criteria:
            print("Missing: filtering criteria")
        usage()
    else:
        bam_filter(infile, outfile, criteria, failed, verbose)
