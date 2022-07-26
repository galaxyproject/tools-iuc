#!/usr/bin/env python

import argparse
import re


AMPLICON_PAT = re.compile(r'.*_(?P<num>\d+).*_(?P<name>L(?:EFT)?|R(?:IGHT)?)')


def write_amplicon_info_file(bed_file, amplicon_info_file):
    amplicon_sets = {}
    for line in bed_file:
        line = line.strip()
        if not line:
            continue
        fields = line.split('\t')
        start = int(fields[1])
        name = fields[3]
        re_match = AMPLICON_PAT.match(name)
        if re_match is None:
            raise ValueError(
                '{} does not match expected amplicon name format'.format(name)
            )
        amplicon_id = int(re_match.group('num'))
        amplicon_set = amplicon_sets.get(amplicon_id, [])
        amplicon_set.append((name, start))
        amplicon_sets[amplicon_id] = amplicon_set

    # write amplicons sorted by number with primers sorted by start position
    for id in sorted(amplicon_sets):
        amplicon_info = '\t'.join(
            [name for name, start in sorted(
                amplicon_sets[id], key=lambda x: x[1]
            )]
        ) + '\n'
        amplicon_info_file.write(amplicon_info)
    amplicon_info_file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Write an amplicon info file for iVar '
                    'from a BED file describing primer positions'
    )
    parser.add_argument(
        'bed_file', type=argparse.FileType(), help='Primer BED file'
    )
    parser.add_argument(
        'amplicon_info_file', type=argparse.FileType('w'),
        help='Output file: amplicon info file in TSV format'
    )
    args = parser.parse_args()

    write_amplicon_info_file(args.bed_file, args.amplicon_info_file)
