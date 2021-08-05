#!/usr/bin/env python

import argparse
import re

AMPLICON_NAME_RE = r'.*_(?P<num>\d+)_[^0-9]*(?P<name>L(?:EFT)?|R(?:IGHT)?)'


def primer_info_to_position(name):
    position = 0
    re_match = re.match(AMPLICON_NAME_RE, name)
    if re_match is None:
        raise ValueError("{} does not match expected amplicon name format".format(name))
    side = re_match.group('name')
    num = re_match.group('num')
    if side == 'RIGHT' or side == 'R':
        position += 1000
    if num is not None:
        position += int(num)
    return position


def write_amplicon_info_file(bed_file, amplicon_info_file):
    amplicon_sets = {}
    amplicon_ids = set()
    for line in bed_file:
        fields = line.strip().split('\t')
        name = fields[3]
        re_match = re.match(AMPLICON_NAME_RE, name)
        if re_match is None:
            raise ValueError("{} does not match expected amplicon name format".format(name))
        amplicon_id = int(re_match.group('num'))
        amplicon_set = amplicon_sets.get(amplicon_id, [])
        amplicon_set.append(name)
        amplicon_sets[amplicon_id] = amplicon_set
        amplicon_ids.add(amplicon_id)

    for id in sorted(list(amplicon_ids)):
        amplicon_info = '\t'.join([name for name in sorted(amplicon_sets[id], key=primer_info_to_position)]) + '\n'
        amplicon_info_file.write(amplicon_info)
    amplicon_info_file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Write an amplicon info file for iVar from a BED file describing primer positions')
    parser.add_argument('bed_file', type=argparse.FileType(), help='Primer BED file')
    parser.add_argument('amplicon_info_file', type=argparse.FileType('w'), help='Output file: amplicon info file in TSV format')
    args = parser.parse_args()

    write_amplicon_info_file(args.bed_file, args.amplicon_info_file)
