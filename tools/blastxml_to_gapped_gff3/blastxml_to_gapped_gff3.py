#!/usr/bin/env python
import argparse
import copy
import logging
import re
import sys

from BCBio import GFF
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(name='blastxml2gff3')

__doc__ = """
BlastXML files, when transformed to GFF3, do not normally show gaps in the
blast hits. This tool aims to fill that "gap".
"""


def blastxml2gff3(blastxml, min_gap=3, trim=False, trim_end=False, include_seq=False):
    from Bio.Blast import NCBIXML
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    blast_records = NCBIXML.parse(blastxml)
    for idx_record, record in enumerate(blast_records):
        # http://www.sequenceontology.org/browser/release_2.4/term/SO:0000343
        match_type = {  # Currently we can only handle BLASTN, BLASTP
            'BLASTN': 'nucleotide_match',
            'BLASTP': 'protein_match',
        }.get(record.application, 'match')

        recid = record.query
        if ' ' in recid:
            recid = recid[0:recid.index(' ')]

        rec = SeqRecord(Seq("ACTG"), id=recid)
        for idx_hit, hit in enumerate(record.alignments):
            for idx_hsp, hsp in enumerate(hit.hsps):
                qualifiers = {
                    "ID": 'b2g.%s.%s.%s' % (idx_record, idx_hit, idx_hsp),
                    "source": "blast",
                    "score": hsp.expect,
                    "accession": hit.accession,
                    "hit_id": hit.hit_id,
                    "length": hit.length,
                    "hit_titles": hit.title.split(' >'),
                }
                if include_seq:
                    qualifiers.update({
                        'blast_qseq': hsp.query,
                        'blast_sseq': hsp.sbjct,
                        'blast_mseq': hsp.match,
                    })

                for prop in ('score', 'bits', 'identities', 'positives',
                             'gaps', 'align_length', 'strand', 'frame',
                             'query_start', 'query_end', 'sbjct_start',
                             'sbjct_end'):
                    qualifiers['blast_' + prop] = getattr(hsp, prop, None)

                desc = hit.title.split(' >')[0]
                qualifiers['description'] = desc[desc.index(' '):]

                # This required a fair bit of sketching out/match to figure out
                # the first time.
                #
                # the match_start location must account for queries and
                # subjecst that start at locations other than 1
                parent_match_start = hsp.query_start - hsp.sbjct_start
                # The end is the start + hit.length because the match itself
                # may be longer than the parent feature, so we use the supplied
                # subject/hit length to calculate the real ending of the target
                # protein.
                parent_match_end = hsp.query_start + hit.length + hsp.query.count('-')

                # If we trim the left end, we need to trim without losing information.
                used_parent_match_start = parent_match_start
                if trim:
                    if parent_match_start < 1:
                        used_parent_match_start = 0

                if trim or trim_end:
                    if parent_match_end > hsp.query_end:
                        parent_match_end = hsp.query_end + 1

                # The ``match`` feature will hold one or more ``match_part``s
                top_feature = SeqFeature(
                    FeatureLocation(used_parent_match_start, parent_match_end),
                    type=match_type, strand=0,
                    qualifiers=qualifiers
                )

                # Unlike the parent feature, ``match_part``s have sources.
                part_qualifiers = {
                    "source": "blast",
                }
                top_feature.sub_features = []
                for idx_part, (start, end, cigar) in \
                        enumerate(generate_parts(hsp.query, hsp.match,
                                                 hsp.sbjct,
                                                 ignore_under=min_gap)):
                    part_qualifiers['Gap'] = cigar
                    part_qualifiers['ID'] = qualifiers['ID'] + ('.%s' % idx_part)

                    # Otherwise, we have to account for the subject start's location
                    match_part_start = parent_match_start + hsp.sbjct_start + start - 1

                    # We used to use hsp.align_length here, but that includes
                    # gaps in the parent sequence
                    #
                    # Furthermore align_length will give calculation errors in weird places
                    # So we just use (end-start) for simplicity
                    match_part_end = match_part_start + (end - start)

                    top_feature.sub_features.append(
                        SeqFeature(
                            FeatureLocation(match_part_start, match_part_end),
                            type="match_part", strand=0,
                            qualifiers=copy.deepcopy(part_qualifiers))
                    )

                rec.features.append(top_feature)
        rec.annotations = {}
        yield rec


def __remove_query_gaps(query, match, subject):
    """remove positions in all three based on gaps in query

    In order to simplify math and calculations...we remove all of the gaps
    based on gap locations in the query sequence::

        Q:ACTG-ACTGACTG
        S:ACTGAAC---CTG

    will become::

        Q:ACTGACTGACTG
        S:ACTGAC---CTG

    which greatly simplifies the process of identifying the correct location
    for a match_part
    """
    prev = 0
    fq = ''
    fm = ''
    fs = ''
    for position in re.finditer('-', query):
        fq += query[prev:position.start()]
        fm += match[prev:position.start()]
        fs += subject[prev:position.start()]
        prev = position.start() + 1
    fq += query[prev:]
    fm += match[prev:]
    fs += subject[prev:]

    return (fq, fm, fs)


def generate_parts(query, match, subject, ignore_under=3):
    region_q = []
    region_m = []
    region_s = []

    (query, match, subject) = __remove_query_gaps(query, match, subject)

    region_start = -1
    region_end = -1
    mismatch_count = 0
    for i, (q, m, s) in enumerate(zip(query, match, subject)):

        # If we have a match
        if m != ' ' or m == '+':
            if region_start == -1:
                region_start = i
                # It's a new region, we need to reset or it's pre-seeded with
                # spaces
                region_q = []
                region_m = []
                region_s = []
            region_end = i
            mismatch_count = 0
        else:
            mismatch_count += 1

        region_q.append(q)
        region_m.append(m)
        region_s.append(s)

        if mismatch_count >= ignore_under and region_start != -1 and region_end != -1:
            region_q = region_q[0:-ignore_under]
            region_m = region_m[0:-ignore_under]
            region_s = region_s[0:-ignore_under]
            yield region_start, region_end + 1, \
                cigar_from_string(region_q, region_m, region_s, strict_m=True)
            region_q = []
            region_m = []
            region_s = []

            region_start = -1
            region_end = -1
            mismatch_count = 0

    yield region_start, region_end + 1, \
        cigar_from_string(region_q, region_m, region_s, strict_m=True)


def _qms_to_matches(query, match, subject, strict_m=True):
    matchline = []

    for (q, m, s) in zip(query, match, subject):
        ret = ''

        if m != ' ' or m == '+':
            ret = '='
        elif m == ' ':
            if q == '-':
                ret = 'D'
            elif s == '-':
                ret = 'I'
            else:
                ret = 'X'
        else:
            log.warn("Bad data: \n\t%s\n\t%s\n\t%s\n" % (query, match, subject))

        if strict_m:
            if ret == '=' or ret == 'X':
                ret = 'M'

        matchline.append(ret)
    return matchline


def _matchline_to_cigar(matchline):
    cigar_line = []
    last_char = matchline[0]
    count = 0
    for char in matchline:
        if char == last_char:
            count += 1
        else:
            cigar_line.append("%s%s" % (last_char, count))
            count = 1
        last_char = char
    cigar_line.append("%s%s" % (last_char, count))
    return ' '.join(cigar_line)


def cigar_from_string(query, match, subject, strict_m=True):
    matchline = _qms_to_matches(query, match, subject, strict_m=strict_m)
    if len(matchline) > 0:
        return _matchline_to_cigar(matchline)
    else:
        return ""


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert Blast XML to gapped GFF3', epilog='')
    parser.add_argument('blastxml', type=argparse.FileType("r"), help='Blast XML Output')
    parser.add_argument('--min_gap', type=int, help='Maximum gap size before generating a new match_part', default=3)
    parser.add_argument('--trim', action='store_true', help='Trim blast hits to be only as long as the parent feature')
    parser.add_argument('--trim_end', action='store_true', help='Cut blast results off at end of gene')
    parser.add_argument('--include_seq', action='store_true', help='Include sequence')
    args = parser.parse_args()

    for rec in blastxml2gff3(**vars(args)):
        if len(rec.features):
            GFF.write([rec], sys.stdout)
