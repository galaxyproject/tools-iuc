import os
import re
import sys

import ngsutils.support
import pysam
try:
    from eta import ETA
except:
    pass


def bam_open(fname, mode='r', *args, **kwargs):
    if fname.lower()[-4:] == '.bam':
        return pysam.Samfile(fname, '%sb' % mode, *args, **kwargs)
    return pysam.Samfile(fname, '%s' % mode, *args, **kwargs)


def bam_pileup_iter(bam, mask=1796, quiet=False, callback=None):
    if not quiet and bam.filename:
        eta = ETA(os.stat(bam.filename).st_size)
    else:
        eta = None

    for pileup in bam.pileup(mask=mask):
        pos = bam.tell()
        bgz_offset = pos >> 16

        if not quiet:
            if callback:
                eta.print_status(bgz_offset, extra=callback(pileup))
            else:
                eta.print_status(bgz_offset, extra='%s:%s' % (bam.getrname(pileup.tid), pileup.pos))

        yield pileup

    if eta:
        eta.done()


def bam_iter(bam, quiet=False, show_ref_pos=False, ref=None, start=None, end=None, callback=None):
    '''
    >>> [x.qname for x in bam_iter(bam_open(os.path.join(os.path.dirname(__file__), 't', 'test.bam')), quiet=True)]
    ['A', 'B', 'E', 'C', 'D', 'F', 'Z']
    '''

    if os.path.exists('%s.bai' % bam.filename):
        # This is an indexed file, so it is ref sorted...
        # Meaning that we should show chrom:pos, instead of read names
        show_ref_pos = True

    eta = None

    if not ref:
        if not quiet and bam.filename:
            eta = ETA(os.stat(bam.filename).st_size)

        for read in bam:
            pos = bam.tell()
            bgz_offset = pos >> 16

            if not quiet and eta:
                if callback:
                    eta.print_status(bgz_offset, extra=callback(read))
                elif (show_ref_pos):
                    if read.tid > -1:
                        eta.print_status(bgz_offset, extra='%s:%s %s' % (bam.getrname(read.tid), read.pos, read.qname))
                    else:
                        eta.print_status(bgz_offset, extra='unmapped %s' % (read.qname))
                else:
                    eta.print_status(bgz_offset, extra='%s' % read.qname)

            yield read

    else:
        working_chrom = None
        if ref in bam.references:
            working_chrom = ref
        elif ref[0:3] == 'chr':
            # compensate for Ensembl vs UCSC ref naming
            if ref[3:] in bam.references:
                working_chrom = ref[3:]

        if not working_chrom:
            raise ValueError('Missing reference: %s' % ref)

        tid = bam.gettid(working_chrom)

        if not start:
            start = 0
        if not end:
            end = bam.lengths[tid]

        if not quiet and bam.filename:
            eta = ETA(end - start)

        for read in bam.fetch(working_chrom, start, end):
            if not quiet and eta:
                if callback:
                    eta.print_status(read.pos - start, extra=callback(read))
                else:
                    eta.print_status(read.pos - start, extra='%s:%s %s' % (bam.getrname(read.tid), read.pos, read.qname))

            yield read

    if eta:
        eta.done()


def bam_batch_reads(bam, quiet=False):
    '''
    Batch mapping for the same reads (qname) together, this way
    they can all be compared/converted together.
    '''
    reads = []
    last = None
    for read in bam_iter(bam, quiet=quiet):
        if last and read.qname != last:
            yield reads
            reads = []
        last = read.qname
        reads.append(read)

    if reads:
        yield reads


bam_cigar = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X']
bam_cigar_op = {
    'M': 0,
    'I': 1,
    'D': 2,
    'N': 3,
    'S': 4,
    'H': 5,
    'P': 6,
    '=': 7,
    'X': 8,
}


def cigar_fromstr(s):
    '''
    >>> cigar_fromstr('10M5I20M')
    [(0, 10), (1, 5), (0, 20)]
    >>> cigar_fromstr('1M1I1D1N1S1H1P')
    [(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1), (6, 1)]
    '''
    ret = []
    spl = re.split('([0-9]+)', s)[1:]
    for op, size in zip(spl[1::2], spl[::2]):
        ret.append((bam_cigar_op[op], int(size)))
    return ret


def cigar_tostr(cigar):
    '''
    >>> cigar_tostr(((0, 10), (1, 5), (0, 20)))
    '10M5I20M'
    >>> cigar_tostr(((0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1), (6, 1)))
    '1M1I1D1N1S1H1P'
    '''

    s = ''

    for op, size in cigar:
        s += '%s%s' % (size, bam_cigar[op])

    return s


def cigar_read_len(cigar):
    '''
    >>> cigar_read_len(cigar_fromstr('8M'))
    8
    >>> cigar_read_len(cigar_fromstr('8M100N8M'))
    16
    >>> cigar_read_len(cigar_fromstr('8M10I8M'))
    26
    >>> cigar_read_len(cigar_fromstr('8M10D8M'))
    16
    '''

    read_pos = 0

    for op, length in cigar:
        if op == 0:  # M
            read_pos += length
        elif op == 1:  # I
            read_pos += length
        elif op == 2:  # D
            pass
        elif op == 3:  # N
            pass
        elif op == 4:  # S
            read_pos += length
        elif op == 5:  # H
            pass
        elif op == 7:  # =
            read_pos += length
        elif op == 8:  # X
            read_pos += length
        else:
            raise ValueError("Unsupported CIGAR operation: %s" % op)

    return read_pos


def read_calc_mismatches(read):
    inserts = 0
    deletions = 0
    indels = 0
    edits = int(read.opt('NM'))
    #
    # NM counts the length of indels
    # We really just care about *if* there is an indel, not the size
    #

    for op, length in read.cigar:
        if op == 1:
            inserts += length
            indels += 1
        elif op == 2:
            deletions += length
            indels += 1

    return edits - inserts - deletions + indels


def _extract_md_matches(md, maxlength):
    md_pos = 0

    while md and md_pos < maxlength:
        # preload a zero so that immediate mismatches will be caught
        # the zero will have no affect otherwise...
        tmp = '0'

        # look for matches
        while md and md[0] in '0123456789':
            tmp += md[0]
            md = md[1:]

        pos = int(tmp)
        if pos > maxlength:
            return (maxlength, '%s%s' % (pos - maxlength, md))
        return (pos, md)


def read_calc_variations(read):
    'see _read_calc_variations'
    for tup in _read_calc_variations(read.pos, read.cigar, read.opt('MD'), read.seq):
        yield tup


def _read_calc_variations(start_pos, cigar, md, seq):
    '''
    For each variation, outputs a tuple: (op, pos, seq)

    op  - operation (0 = mismatch, 1 = insert, 2 = deletion) (like CIGAR)
    pos - 0-based position of the variation (relative to reference)
    seq - the base (or bases) involved in the variation
          for mismatch or insert, this is the sequence inserted
          for deletions, this is the reference sequence that was removed

    MD is the mismatch string. Not all aligners include the tag. If your aligner
    doesn't include this, then you'll need to add it, or use a different function
    (see: read_calc_mismatches_gen).

    Special care must be used to handle RNAseq reads that cross
    an exon-exon junction.

    Also: MD is a *really* dumb format that can't be read correctly with
          a regex. It must be processed in concert with the CIGAR alignment
          in order to catch all edge cases. Some implementations insert 0's
          at the end of inserts / deltions / variations to make parsing easier
          but not everyone follows this. Look at the complex examples: the
          CIGAR alignment may show an insert, but the MD just shows all matches.

    Examples: See: http://davetang.org/muse/2011/01/28/perl-and-sam/
              Also from CCBB actual mappings and manual altered (shortened,
              made more complex)
              (doctests included)

    Match/mismatch
    CIGAR: 36M
    MD:Z:  1A0C0C0C1T0C0T27
    MD:Z:  1ACCC1TCT27 (alternative)
                   1         2
          123456789012345678901234567890123456
    ref:  CGATACGGGGACATCCGGCCTGCTCCTTCTCACATG
           XXXX XXX
    read: CACCCCTCTGACATCCGGCCTGCTCCTTCTCACATG
          MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
          -ACCC-TCT---------------------------
    >>> list(_read_calc_variations(1, [(0,36)], '1A0C0C0C1T0C0T27', 'CACCCCTCTGACATCCGGCCTGCTCCTTCTCACATG'))
    [(0, 2, 'A'), (0, 3, 'C'), (0, 4, 'C'), (0, 5, 'C'), (0, 7, 'T'), (0, 8, 'C'), (0, 9, 'T')]

    Insert
    CIGAR: 6M1I29M
    MD:Z: 0C1C0C1C0T0C27
          C1CC1CTC27 (alt)
                    1         2
          123456^789012345678901234567890123456
    ref:  CACCCC^TCTGACATCCGGCCTGCTCCTTCTCACAT
          X XX X|XX
    read: GAGACGGGGTGACATCCGGCCTGCTCCTTCTCACAT
          MMMMMMIMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
          G-GA-GGGG---------------------------
    >>> list(_read_calc_variations(1, [(0,6), (1,1), (0, 29)], '0C1C0C1C0T0C27', 'GAGACGGGGTGACATCCGGCCTGCTCCTTCTCACAT'))
    [(0, 1, 'G'), (0, 3, 'G'), (0, 4, 'A'), (0, 6, 'G'), (1, 7, 'G'), (0, 7, 'G'), (0, 8, 'G')]
    >>> list(_read_calc_variations(1, [(0,6), (1,1), (0, 29)], 'C1CC1CTC27', 'GAGACGGGGTGACATCCGGCCTGCTCCTTCTCACAT'))
    [(0, 1, 'G'), (0, 3, 'G'), (0, 4, 'A'), (0, 6, 'G'), (1, 7, 'G'), (0, 7, 'G'), (0, 8, 'G')]


    Deletion
    CIGAR: 9M9D27M
    MD:Z: 2G0A5^ATGATGTCA27
          2GA5^ATGATGTCA27 (alt)
    ref:  AGGAATGGGATGATGTCAGGGGTTCCAGGTGGAGACGAGGACTCC
            XX     ^^^^^^^^^
    read: AGTGATGGG^^^^^^^^^GGGGTTCCAGGTGGAGACGAGGACTCC
          MMMMMMMMMDDDDDDDDDMMMMMMMMMMMMMMMMMMMMMMMMMMM
          --TG-----ATGATGTCA---------------------------
    >>> list(_read_calc_variations(1, [(0,9), (2,9), (0, 27)], '2G0A5^ATGATGTCA27', 'AGTGATGGGGGGGTTCCAGGTGGAGACGAGGACTCC'))
    [(0, 3, 'T'), (0, 4, 'G'), (2, 10, 'ATGATGTCA')]


    Complex
    CIGAR: 9M9D11M1I15M
    MD:Z: 2G0A5^ATGATGTCAA26
    MD:Z: 2G0A5^ATGATGTCA0G26 (alt)
                   1         2         3         4
    pos:  123456789012345678901234567890123456789012345
    ref:  AGGAATGGGATGATGTCAGGGGTTCCAGG^GGAGACGAGGACTCC
            XX     ^^^^^^^^^X          |
    read: AGTGATGGG^^^^^^^^^AGGGTTCCAGGTGGAGACGAGGACTCC
          MMMMMMMMMDDDDDDDDDMMMMMMMMMMMMMMMMMMMMMMMMMMM
          --TG-----ATGATGTCAG----------T---------------
    >>> list(_read_calc_variations(1, [(0,9), (2,9), (0,11), (1,1), (0,15)], '2G0A5^ATGATGTCAA26', 'AGTGATGGGGGGGTTCCAGGTGGAGACGAGGACTCC'))
    [(0, 3, 'T'), (0, 4, 'G'), (2, 10, 'ATGATGTCA'), (0, 19, 'G'), (1, 30, 'T')]


    Complex example - inserts aren't separately handled by MD, only visible in CIGAR
    CIGAR: 14M2D16M3I42M
    MD:Z:  14^TC58
                   1         2         3            4         5         6         7
    pos:  12345678901234567890123456789012^^^345678901234567890123456789012345678901234567
    ref:  caagtatcaccatgtcaggcatttttttcatt^^^tttgtagagagagaagacttgctatgttgcccaagctggcct
                        ^^                |||
    read: CAAGTATCACCATG^^AGGCATTTTTTTCATTTGGTTTGTAGAGAGAGAAGACTTGCTATGTTGCCCAAGCTGGCCT
          MMMMMMMMMMMMMMDDMMMMMMMMMMMMMMMMIIIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
          --------------tc----------------TGG------------------------------------------
    >>> list(_read_calc_variations(1, [(0,14), (2,2), (0,16), (1,3), (0,42)], '14^TC58', 'CAAGTATCACCATGAGGCATTTTTTTCATTTGGTTTGTAGAGAGAGAAGACTTGCTATGTTGCCCAAGCTGGCCT'))
    [(2, 15, 'TC'), (1, 33, 'TGG')]


    Complex example 2:
    CIGAR: 41M3I10M1I5M1I2M2I10M
    MD:Z:  44C2C6T6T6
                   1         2         3         4            5             6
    pos:  12345678901234567890123456789012345678901^^^2345678901^23456^78^^9012345678
    ref:  AGGGTGGCGAGATCGATGACGGCATTGGCGATGGTGATCTT^^^GAGCCACATG^CGGTC^GC^^GGATCTCCAG
                                                   |||   X  X   |   X |  ||   X
    read: AGGGTGGCGAGATCGATGACGGCATTGGCGATGGTGATCTTTTAGAGACATATGCCGGACGGCGTGGAGCTCCAG
          MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMIIIMMMMMMMMMMIMMMMMIMMIIMMMMMMMMMM
          -----------------------------------------tta---A--T---c---A----gt---G------


    13M 28M 3I 10M 1I 5M 1I 2M 2I 10M
    >>> list(_read_calc_variations(1, [(0, 41), (1, 3), (0, 10), (1, 1), (0, 5), (1, 1), (0, 2), (1, 2), (0, 10)], '44C2C6T6T6', 'AGGGTGGCGAGATCGATGACGGCATTGGCGATGGTGATCTTTTAGAGACATATGCCGGACGGCGTGGAGCTCCAG'))
    [(1, 42, 'TTA'), (0, 45, 'A'), (0, 48, 'T'), (1, 52, 'C'), (0, 55, 'A'), (1, 57, 'G'), (1, 59, 'GT'), (0, 62, 'G')]


    Splice junction example:
    CIGAR: 62M100N13M
    MD:Z:  2T27C44
                                                                                 1      1
                   1         2         3         4         5         6           6      7
    pos:  12345678901234567890123456789012345678901234567890123456789012| [100] |3456789012345
    ref:  CCTCATGACCAGCTTGTTGAAGAGATCCGACATCAAGTGCCCACCTTGGCTCGTGGCTCTCA|-------|CTTGCTCCTGCTC
            X                           X
    read: CCGCATGACCAGCTTGTTGAAGAGATCCGATATCAAGTGCCCACCTTGGCTCGTGGCTCTCA|-------|CTTGCTCCTGCTC
          --G---------------------------T-----------------------------------------------------

    >>> list(_read_calc_variations(1, [(0,62), (4,100), (0,13)], '2T27C44', 'CCGCATGACCAGCTTGTTGAAGAGATCCGATATCAAGTGCCCACCTTGGCTCGTGGCTCTCACTTGCTCCTGCTC'))
    [(0, 3, 'G'), (0, 31, 'T')]


    Splice junction example 2:
    CIGAR: 13M100N28M3I10M1I5M1I2M2I10M
    MD:Z:  44C2C6T6T6
                                      1         1         1            1             1
                   1                  2         3         4            5             6
    pos:  1234567890123| [100] |4567890123456789012345678901^^^2345678901^23456^78^^9012345678
    ref:  AGGGTGGCGAGAT|-------|CGATGACGGCATTGGCGATGGTGATCTT^^^GAGCCACATG^CGGTC^GC^^GGATCTCCAG
                                                            |||   X  X   |   X |  ||   X
    read: AGGGTGGCGAGAT|-------|CGATGACGGCATTGGCGATGGTGATCTTTTAGAGACATATGCCGGACGGCGTGGAGCTCCAG
          MMMMMMMMMMMMM         MMMMMMMMMMMMMMMMMMMMMMMMMMMMIIIMMMMMMMMMMIMMMMMIMMIIMMMMMMMMMM
          -------------         ----------------------------tta---A--T---c---A----gt---G------

    13M 100N 28M 3I 10M 1I 5M 1I 2M 2I 10M
    >>> list(_read_calc_variations(1, [(0, 13), (3, 100), (0, 28), (1, 3), (0, 10), (1, 1), (0, 5), (1, 1), (0, 2), (1, 2), (0, 10)], '44C2C6T6T6', 'AGGGTGGCGAGATCGATGACGGCATTGGCGATGGTGATCTTTTAGAGACATATGCCGGACGGCGTGGAGCTCCAG'))
    [(1, 142, 'TTA'), (0, 145, 'A'), (0, 148, 'T'), (1, 152, 'C'), (0, 155, 'A'), (1, 157, 'G'), (1, 159, 'GT'), (0, 162, 'G')]


    Splice junction example 2A:
    CIGAR: 13M100N7M2D19M3I10M1I5M1I2M2I10M
    MD:Z:  9A10^GG22C2C6T6T6
                                      1         1         1            1             1
                   1                  2         3         4            5             6
    pos:  1234567890123| [100] |4567890123456789012345678901^^^2345678901^23456^78^^9012345678
    ref:  AGGGTGGCGAGAT|-------|CGATGACGGCATTGGCGATGGTGATCTT^^^GAGCCACATG^CGGTC^GC^^GGATCTCCAG
                                       ^^                   |||   X  X   |   X |  ||   X
    read: AGGGTGGCGCGAT|-------|CGATGAC^^CATTGGCGATGGTGATCTTTTAGAGACATATGCCGGACGGCGTGGAGCTCCAG
          MMMMMMMMMMMMM         MMMMMMMDDMMMMMMMMMMMMMMMMMMMIIIMMMMMMMMMMIMMMMMIMMIIMMMMMMMMMM
          ---------C---         ----------------------------tta---A--T---c---A----gt---G------
          .........A...         .......GG...................   ...C..C... ...T. ..  ...T......
              9    A        10        ^GG             22          C 2C   6   T     6   T   6

    >>> list(_read_calc_variations(1, [(0, 13), (3, 100), (0, 7), (2, 2), (0, 19), (1, 3), (0, 10), (1, 1), (0, 5), (1, 1), (0, 2), (1, 2), (0, 10)], '9A10^GG22C2C6T6T6', 'AGGGTGGCGCGATCGATGACCATTGGCGATGGTGATCTTTTAGAGACATATGCCGGACGGCGTGGAGCTCCAG'))
    [(0, 10, 'C'), (2, 121, 'GG'), (1, 142, 'TTA'), (0, 145, 'A'), (0, 148, 'T'), (1, 152, 'C'), (0, 155, 'A'), (1, 157, 'G'), (1, 159, 'GT'), (0, 162, 'G')]

    Real Example
    242_1071_1799_B1
    CIGAR: 42M10I3M1D9M1D11M
    MD:Z:  27G16A0^T6C2^T1C9
                   1         2         3         4                   5         6         7
    pos:  123456789012345678901234567890123456789012          345678901234567890123456789012345
    ref:  ACTGAGAAACCCAACCCTCTGAGACCAGCACACCCCTTTCAA^^^^^^^^^^GCATGTTCCTCCCTCCCCTTCTTTG
                                     X                          X^      X  ^ X
    read: ACTGAGAAACCCAACCCTCTGAGACCAACACACCCCTTTCAACACATTTTTGGCC^GTTCCTGCC^CGCCTTCTTTG
          MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMIIIIIIIIIIMMMDMMMMMMMMMDMMMMMMMMMMM
          ---------------------------A--------------^^^^^^^^^^--CT------G--T-G---------

    >>> list(_read_calc_variations(1, [(0,42), (1,10), (0, 3), (2, 1), (0, 9), (2, 1), (0, 11)], '27G16A0^T6C2^T1C9', 'ACTGAGAAACCCAACCCTCTGAGACCAACACACCCCTTTCAACACATTTTTGGCCGTTCCTGCCCGCCTTCTTTG',  ))
    [(0, 28, 'A'), (1, 43, 'CACATTTTTG'), (0, 45, 'C'), (2, 46, 'T'), (0, 53, 'G'), (2, 56, 'T'), (0, 58, 'G')]


    Real example 2
    577_1692_891_A1
    CIGAR: 34M100N39M2I
    MD:Z:  3T69
                                                          1         1         1         1
                   1         2         3                  4         5         6         7
    pos:  1234567890123456789012345678901234| [100] |567890123456789012345678901234567890123
    ref:  GGATTCTTCCCACTGGGTCGATGTTGTTTGTGAT|-------|CTGAGAGAGAGTTGCATCTGCACATGCTTTCCTGGCGTC^^

    read: GGAATCTTCCCACTGGGTCGATGTTGTTTGTGAT|-------|CTGAGAGAGAGTTGCATCTGCACATGCTTTCCTGGCGTCTC
          MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM  NNNNN  MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMII
          ---A------------------------------         ---------------------------------------TC

    >>> list(_read_calc_variations(1, [(0,34), (3,100), (0, 39), (1, 2)], '3T69', 'GGAATCTTCCCACTGGGTCGATGTTGTTTGTGATCTGAGAGAGAGTTGCATCTGCACATGCTTTCCTGGCGTCTC',  ))
    [(0, 4, 'A'), (1, 174, 'TC')]

    '''

    ref_pos = start_pos
    read_pos = 0

    for op, length in cigar:
        if md and md[0] == '0':
            md = md[1:]
        # sys.stderr.write('%s, %s, %s\n' %(op, length, md))
        if op == 0:  # M
            # how far in the chunk are we? (do *not* update ref_pos until end)
            md_pos = 0
            last = None
            while md and md_pos < length:
                if last == (op, length, md):
                    sys.stderr.write('\nInfinite loop in variant finding!\nPos: %s\nCIGAR: (%s, %s)\n' % (ref_pos, op, length))
                    sys.exit(1)
                last = (op, length, md)
                # sys.stderr.write('%s, %s, %s\n' %(op, length, md))
                chunk_size, md = _extract_md_matches(md, length - md_pos)
                # sys.stderr.write('   -> %s, %s\n' %(chunk_size, md))
                md_pos += chunk_size

                # look for mismatches
                while md_pos < length and md and md[0] not in '0123456789^':
                    yield (op, ref_pos + md_pos, seq[read_pos + md_pos])
                    md = md[1:]

                    md_pos += 1

            ref_pos += length
            read_pos += length

        elif op == 1:  # I
            # nothing in MD about inserts...
            yield (op, ref_pos, seq[read_pos:read_pos + length])
            read_pos += length

        elif op == 2:  # D
            # prefixed with '^' and includes all of the removed bases
            if md[0] == '^':
                md = md[1:]
            yield (op, ref_pos, md[:length])
            md = md[length:]
            ref_pos += length

        elif op == 3:  # N
            ref_pos += length


def read_calc_mismatches_gen(ref, read, chrom):
    start = read.pos
    ref_pos = 0
    read_pos = 0

    for op, length in read.cigar:
        if op == 1:
            yield ref_pos, op, None
            read_pos += length
        elif op == 2:
            yield ref_pos, op, None
            ref_pos += length
        elif op == 3:
            ref_pos += length
        elif op == 0:
            refseq = ref.fetch(chrom, start + ref_pos, start + ref_pos + length)
            if not refseq:
                raise ValueError("Reference '%s' not found in FASTA file: %s" % (chrom, ref.filename))
            cur_pos = start + ref_pos
            for refbase, readbase in zip(refseq.upper(), read.seq[read_pos:read_pos + length].upper()):
                if refbase != readbase:
                    yield op, cur_pos, readbase
                cur_pos += 1
            ref_pos += length
            read_pos += length
        else:
            raise ValueError("Unsupported CIGAR operation: %s" % op)


def read_calc_mismatches_ref(ref, read, chrom):
    edits = 0

    for op, pos, base in read_calc_mismatches_gen(ref, read, chrom):
        edits += 1

    return edits


__region_cache = {}


def region_pos_to_genomic_pos(name, start, cigar):
    '''
        converts a junction position to a genomic location given a junction
        ref name, the junction position, and the cigar alignment.

        returns: (genomic ref, genomic pos, genomic cigar)

    >>> region_pos_to_genomic_pos('chr1:1000-1050,2000-2050,3000-4000', 25, [(0, 100)])
    ('chr1', 1025, [(0, 25), (3, 950), (0, 50), (3, 950), (0, 25)])

    >>> region_pos_to_genomic_pos('chr1:1000-1050,1050-1200', 25, [(0, 100)])
    ('chr1', 1025, [(0, 25), (3, 0), (0, 75)])

    >>> region_pos_to_genomic_pos('chr3R:17630851-17630897,17634338-17634384', 17, [(0, 39)])
    ('chr3R', 17630868, [(0, 29), (3, 3441), (0, 10)])

    >>> region_pos_to_genomic_pos('chr1:1000-1050,2000-2050', 25, [(4, 25), (0, 50), (4, 25)])
    ('chr1', 1025, [(4, 25), (0, 25), (3, 950), (0, 25), (4, 25)])

    >>> region_pos_to_genomic_pos('chr1:1000-1050,2000-2050', 25, [(5, 25), (0, 75)])
    ('chr1', 1025, [(5, 25), (0, 25), (3, 950), (0, 50)])

    >>> region_pos_to_genomic_pos('chr1:1000-1050,2000-2050', 25, [(5, 25), (0, 75)])
    ('chr1', 1025, [(5, 25), (0, 25), (3, 950), (0, 50)])

    >>> region_pos_to_genomic_pos('chr7:16829153-16829246,16829246-16829339', 62, cigar_fromstr('83M18S'))
    ('chr7', 16829215, [(0, 31), (3, 0), (0, 52), (4, 18)])


    '''

    if name in __region_cache:
        chrom, fragments = __region_cache[name]
    else:
        c1 = name.split(':')
        chrom = c1[0]

        fragments = []
        for fragment in c1[1].split(','):
            s, e = fragment.split('-')
            fragments.append((int(s), int(e)))

        __region_cache[name] = (chrom, fragments)

    chr_cigar = []
    chr_start = fragments[0][0]

    read_start = int(start)

    frag_idx = 0
    frag_start = 0
    frag_end = 0

    for i, (s, e) in enumerate(fragments):
        if chr_start + read_start < e:
            chr_start += read_start
            frag_idx = i
            frag_start = s
            frag_end = e
            break

        else:
            chr_start += (e - s)
            read_start -= (e - s)

    cur_pos = chr_start

    for op, length in cigar:
        if op in [1, 4, 5]:
            chr_cigar.append((op, length))

        elif op in [0, 2]:
            if cur_pos + length <= frag_end:
                cur_pos += length
                chr_cigar.append((op, length))

            else:
                while cur_pos + length > frag_end:
                    if frag_end - cur_pos > 0:
                        chr_cigar.append((op, frag_end - cur_pos))
                        length -= (frag_end - cur_pos)
                    cur_pos = frag_end
                    frag_idx += 1
                    if len(fragments) <= frag_idx:
                        print 'ERROR converting: ', name, fragments
                        return (chrom, 0, chr_cigar)
                    frag_start, frag_end = fragments[frag_idx]
                    chr_cigar.append((3, frag_start - cur_pos))
                    cur_pos = frag_start

                cur_pos = cur_pos + length
                chr_cigar.append((op, length))
        else:
            print "Unsupported CIGAR operation (%s)" % bam_cigar[op]
            sys.exit(1)

    return (chrom, chr_start, chr_cigar)


def is_junction_valid(cigar, min_overlap=4):
    '''
        Does the genomic cigar alignment represent a 'good' alignment.
        Used for checking junction->genome alignments

        1) the alignment must not start at a splice junction
        2) the alignment must not start or end with an overhang
        3) the alignment must overhang the splice junction by min_overlap (4)

        |     Exon1       |     Intron     |      Exon2       |
        |-----------------|oooooooooooooooo|------------------|
                                            XXXXXXXXXXXXXXXXXXXXXXXX (bad 1)
      XXXXXXXXXXXXX (bad 2)                           XXXXXXXXXXXXXXXX (bad 2)
                        XX-----------------XXXXXXXXXXXXXXXXX (bad 3)

    >>> is_junction_valid(cigar_fromstr('1000N40M'))
    (False, 'Starts at gap (1000N40M)')

    >>> is_junction_valid(cigar_fromstr('100M'))
    (False, "Doesn't cover junction")

    >>> is_junction_valid(cigar_fromstr('100M1000N3M'), 4)
    (False, "Too short overlap at 3' (100M1000N3M)")

    >>> is_junction_valid(cigar_fromstr('2M1000N100M'), 4)
    (False, "Too short overlap at 5' (2M1000N100M)")

    >>> is_junction_valid(cigar_fromstr('4M1000N100M'), 4)
    (True, '')

    >>> is_junction_valid(cigar_fromstr('100M1000N4M'), 4)
    (True, '')

    >>> is_junction_valid(cigar_fromstr('4M0N100M'), 4)
    (True, '')

    '''
    first = True
    pre_gap = True

    pre_gap_count = 0
    post_gap_count = 0

    has_gap = False

    for op, length in cigar:
        # mapping can't start at a gap
        if first and op == 3:
            return (False, 'Starts at gap (%s)' % cigar_tostr(cigar))
        first = False

        if op == 3:
            pre_gap = False
            post_gap_count = 0
            has_gap = True

        elif pre_gap:
            pre_gap_count += length
        else:
            post_gap_count += length

        # mapping must start with more than min_overlap base match

    if not has_gap:
        return (False, "Doesn't cover junction")
    elif pre_gap_count < min_overlap:
        return (False, "Too short overlap at 5' (%s)" % cigar_tostr(cigar))
    elif post_gap_count < min_overlap:
        return (False, "Too short overlap at 3' (%s)" % cigar_tostr(cigar))

    return True, ''


def read_alignment_fragments_gen(read):
    '''
    Takes a read and returns the start and end positions for each match/mismatch region.
    This will let us know where each read alignment "touches" the genome.
    '''

    for start, end in _read_alignment_fragments_gen(read.pos, read.cigar):
        yield (start, end)


def _read_alignment_fragments_gen(pos, cigar):
    ''' Test-able version of read_alignment_fragments_gen
    >>> list(_read_alignment_fragments_gen(1, cigar_fromstr('50M')))
    [(1, 51)]

    >>> list(_read_alignment_fragments_gen(1, cigar_fromstr('25M100N25M')))
    [(1, 26), (126, 151)]

    >>> list(_read_alignment_fragments_gen(1, cigar_fromstr('20M1D4M100N10M5I10M')))
    [(1, 21), (22, 26), (126, 136), (136, 146)]
    '''
    ref_pos = pos

    for op, length in cigar:
        if op == 0:
            yield (ref_pos, ref_pos + length)
            ref_pos += length
        elif op == 1:
            pass
        elif op == 2:
            ref_pos += length
        elif op == 3:
            ref_pos += length
        else:
            raise ValueError("Unsupported CIGAR operation: %s" % op)


def read_cigar_at_pos(cigar, qpos, is_del):
    '''
    Returns the CIGAR operation for a given read position

    qpos is the 0-based index of the base in the read
    '''
    pos = 0
    returnnext = False
    for op, length in cigar:
        if returnnext:
            return op
        if op == 0:
            pos += length
        elif op == 1:
            pos += length
        elif op == 2:
            pass
        elif op == 3:
            pass
        elif op == 4:
            pos += length
        elif op == 5:
            pass
        else:
            raise ValueError("Unsupported CIGAR operation: %s" % op)

        if pos > qpos:
            return op
        if is_del and pos == qpos:
            returnnext = True

    return None


def cleancigar(cigar):
    '''
    Cleans a CIGAR alignment to remove zero length operations

    >>> cigar_tostr(cleancigar(cigar_fromstr('31M0N52M18S')))
    '83M18S'

    '''
    newcigar = []

    changed = False
    zero = False

    for op, size in cigar:
        if size > 0:
            if zero and newcigar:
                if newcigar[-1][0] == op:
                    newsize = newcigar[-1][1] + size
                    newcigar[-1] = (op, newsize)
                zero = 0
            else:
                newcigar.append((op, size))
        else:
            changed = True
            zero = True

    if changed:
        return newcigar

    return None


def read_cleancigar(read):
    '''
    Replaces the CIGAR string for a read to remove any operations that are zero length.

    This happens to the read in-place
    '''
    if read.is_unmapped:
        return False

    newcigar = []

    newcigar = cleancigar(read.cigar)

    if newcigar:
        read.cigar = newcigar
        return True

    return False


def read_to_unmapped(read, ref=None):
    '''
    Converts a read from mapped to unmapped.

    Sets the 'ZR' tag to indicate the original ref/pos/cigar (if ref is passed)
    '''

    newread = pysam.AlignedRead()

    if ref:
        tags = [('ZR', '%s:%s:%s' % (ref, read.pos, cigar_tostr(read.cigar)))]

    newread.is_unmapped = True
    newread.mapq = 0
    newread.tlen = 0
    newread.pos = -1
    newread.pnext = -1
    newread.rnext = -1
    newread.tid = -1

    newread.qname = read.qname

    if read.is_paired:
        newread.is_paired = True

    if not read.is_unmapped and read.is_reverse:
        newread.seq = ngsutils.support.revcomp(read.seq)
        newread.qual = read.qual[::-1]
    else:
        newread.seq = read.seq
        newread.qual = read.qual

    newread.tags = tags

    return newread


if __name__ == '__main__':
    import doctest
    doctest.testmod()
