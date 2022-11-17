#!/usr/bin/env python

# This code is very loosely based on
# the FusionCatcher sam2psl script here:
#  https://github.com/ndaniel/fusioncatcher/blob/master/bin/sam2psl.py

import argparse
import sys


def parse_cigar(c):
    r = []
    d = ''
    mismatches_x = 0
    c = c.upper()
    for a in c:
        if a.isdigit():
            d = d + a
        elif a in ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X']:
            dd = int(d)
            r.append((a, dd))
            if a == 'X':
                mismatches_x = mismatches_x + dd
            d = ''
        else:
            msg = "Error: unknown CIGAR: %s\n" % str(c)
            sys.exit(msg)
    return (r, mismatches_x)


def blocks(cigar, ig=0):
    # Returns block of matches.  Hard clipping is
    # converted to soft clipping index on read.
    ir = 0
    rr = []
    rg = []
    match = 0
    mismatch = 0
    mismatch_x = 0
    mismatch_clip = 0
    insert_query = 0
    insert_query_count = 0
    insert_ref = 0
    insert_ref_count = 0
    # Sum of lengths of the M/I/S/=/X operations
    # will equal the length of SEQ.
    seq_len = 0
    (cig, mismatch_x) = parse_cigar(cigar)
    mismatch = mismatch_x
    for e in cig:
        if e[0] in ('S', 'H'):
            ir = ir + e[1]
            mismatch_clip = mismatch_clip + e[1]
            seq_len = seq_len + e[1]
        elif e[0] in ('I',):
            ir = ir + e[1]
            mismatch = mismatch + e[1]
            insert_query = insert_query + e[1]
            insert_query_count = insert_query_count + 1
            seq_len = seq_len + e[1]
        elif e[0] in ('X'):
            ir = ir + e[1]
            ig = ig + e[1]
            mismatch = mismatch + e[1]
            mismatch_x = mismatch_x + e[1]
            seq_len = seq_len + e[1]
        elif e[0] in ('M', '='):
            rr.append((ir, ir + e[1]))
            rg.append((ig, ig + e[1]))
            ir = ir + e[1]
            ig = ig + e[1]
            match = match + e[1]
            seq_len = seq_len + e[1]
        elif e[0] in ('D', 'N', 'P'):
            ig = ig + e[1]
            insert_ref = insert_ref + e[1]
            insert_ref_count = insert_ref_count + 1
    return (rr, rg, match, mismatch, mismatch_clip, mismatch_x, insert_ref, insert_ref_count, insert_query, insert_query_count, seq_len)


def get_psl(sam, lens):
    # Returns PSL coordinates.
    if sam and sam[1].isdigit():
        unmapped = True if int(sam[1]) & 0x4 else False
        if (not unmapped) and sam[2] != '*' and sam[5] != '*' and sam[0] != '*':
            # Initialize psl_items to those
            # that constitute a PSL empty line.
            psl_items = ['0', '0', '0', '0', '0', '0', '0', '0', '+', 's', '0', '0', '0', 'r', '0', '0', '0', '0', ',', ',', ',']
            # Read sequence length.
            psl_items[14] = lens.get(sam[2], 0)
            # Reference name.
            psl_items[13] = sam[2]
            # Read name.
            psl_items[9] = sam[0]
            # Strand.
            psl_items[8] = "-" if int(sam[1]) & 0x10 else '+'
            # Start position.
            psl_items[15] = int(sam[3]) - 1
            (interval_query, interval_ref, match, mismatch, mismatch_clip, mismatch_x, insert_ref, insert_ref_count, insert_query, insert_query_count, seq_len) = blocks(sam[5], ig=psl_items[15])
            # Read sequence length.
            if sam[9] != '*' and sam[5].find('H') == -1:
                psl_items[10] = len(sam[9])
            else:
                # The length of SEQ is the sum of the
                # lengths of the M/I/S/=/X operations.
                psl_items[10] = seq_len
            psl_items[4] = insert_query_count
            psl_items[5] = insert_query
            psl_items[6] = insert_ref_count
            psl_items[7] = insert_ref
            # Extract the mismatches using tag
            # NM:i (NM is mismatches per reads).
            tag_nm_i = [e.partition("NM:i:")[2] for e in sam[11:] if e.startswith('NM:i:')]
            if not tag_nm_i:
                # NM is not ideal because it is mismatches
                # per # fragment and not per read, but is
                # etter than nothing.
                tag_nm_i = [e.partition("nM:i:")[2] for e in sam[11:] if e.startswith('nM:i:')]
            tag_nm_i = int(tag_nm_i[0]) if tag_nm_i else 0
            if tag_nm_i > float(0.90) * seq_len:
                tag_nm_i = 0
            # Compute the matches and mismatches (include the
            # clipping as mismatches).
            psl_items[0] = match
            psl_items[1] = mismatch
            if interval_query:
                psl_items[11] = interval_query[0][0]
                psl_items[12] = interval_query[-1][1]
                psl_items[16] = interval_ref[-1][1]
                psl_items[17] = len(interval_query)
                # BLAT always gives the coordinates as
                # everything is mapped on the forwward
                # strand.
                psl_items[18] = ','.join([str(e[1] - e[0]) for e in interval_query]) + ','
                psl_items[19] = ','.join([str(e[0]) for e in interval_query]) + ','
                psl_items[20] = ','.join([str(e[0]) for e in interval_ref]) + ','
            return map(str, psl_items)
    else:
        return None


def to_psl(input_file, psl_file):
    # Convert a SAM file to PSL format.
    header = dict()
    with open(input_file, 'r') as fin, open(psl_file, 'w') as fou:
        for i, sam_line in enumerate(fin):
            sam_line = sam_line.rstrip('\r\n')
            if sam_line:
                sam_items = sam_line.split('\t')
                if sam_items[0].startswith('@'):
                    if sam_items[0].startswith('@SQ') and sam_items[1].startswith('SN:') and sam_items[2].startswith('LN:'):
                        k = sam_items[1][3:]
                        v = int(sam_items[2][3:])
                        header[k] = v
                else:
                    psl_items = get_psl(sam_items, header)
                    if psl_items is not None:
                        fou.write('%s\n' % '\t'.join(psl_items))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("--input_file", action="store", dest="input_file", help="Input file in SAM format.")
    parser.add_argument("--output_file", action="store", dest="output_file", help="Output file in PSL format.")

    args = parser.parse_args()

    to_psl(args.input_file, args.output_file)
