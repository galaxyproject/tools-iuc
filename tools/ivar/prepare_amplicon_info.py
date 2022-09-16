#!/usr/bin/env python

# extends ivar trim's amplicon info parsing abilities
# to include calculation of amplicon regions from
# sets of nested (more than two) primers

import sys


# parse primers and their start positions from BED file
primer_starts = {}
with open(sys.argv[1]) as i:
    for line in i:
        line = line.strip()
        if not line:
            continue
        f = line.split('\t')
        try:
            if f[5] == '+':
                primer_starts[f[3]] = int(f[1])
            elif f[5] == '-':
                primer_starts[f[3]] = int(f[2]) - 1
            else:
                raise ValueError()
        except (IndexError, ValueError):
            sys.exit(
                'Primer BED file needs to be TAB-separated with the '
                'following columns: '
                'chrom, chromStart, chromEnd, name, score, strand, '
                'where "chromStart", "chromEnd" need to be integer values '
                'and "strand" needs to be either "+" or "-".'
            )

# parse amplicon info and record outer primer names
with open(sys.argv[2]) as i:
    ret_lines = []
    for line in i:
        line = line.strip()
        if not line:
            continue
        first = last = None
        for pname in line.split('\t'):
            try:
                primer_start = primer_starts[pname]
            except KeyError:
                sys.exit(
                    'Amplicon info with primer name not found in '
                    f'primer BED file: "{pname}"'
                )
            if first is None or primer_start < primer_starts[first]:
                first = pname
            if last is None or primer_start > primer_starts[last]:
                last = pname
        if first == last:
            sys.exit(
                line
                + 'is not a proper amplicon info line.'
            )
        ret_lines.append(f'{first}\t{last}\n')

# write amended amplicon info
with open(sys.argv[3], 'w') as o:
    o.writelines(ret_lines)
