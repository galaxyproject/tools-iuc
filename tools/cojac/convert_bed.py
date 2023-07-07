#!/usr/bin/env python

import argparse
import re

class Scheme:
    amplicon_pat = re.compile(
        r'(?P<prefix>(.*_)*)(?P<num>\d+).*_(?P<name>L(?:EFT)?|R(?:IGHT)?)'
    )

    @classmethod
    def infer_from_primer_scheme(cls, primer_scheme):
        amplicons = {}
        prefixes_seen = set()

        for primer_dat, ref_id, pool_id in cls.read_primer_bed(primer_scheme):
            name = primer_dat[2]
            re_match = cls.amplicon_pat.match(name)
            if re_match is None:
                raise ValueError(
                    '{} does not match expected amplicon name format'.format(name)
                )
            prefix = re_match.group('prefix')[:-1]
            prefixes_seen.add(prefix)
            amplicon_id = int(re_match.group('num'))
            if amplicon_id in amplicons:
                amplicon = amplicons[amplicon_id]
                if amplicon[1] != ref_id:
                    raise ValueError('ref error')
                if amplicon[2] != pool_id:
                    raise ValueError('pool_error')
                amplicon[0].append(primer_dat)
            else:
                amplicons[amplicon_id] = ([primer_dat], ref_id, pool_id)

        if len(prefixes_seen) == 1 and prefix:
            scheme_name = prefix
        else:
            scheme_name = None

        return cls(amplicons, scheme_name)

    @classmethod
    def from_primers_and_amplicons(
        cls, primer_scheme, amplicon_info, scheme_name=None
    ):
        primer_amplicon_mapping = {}
        amplicon_id = 1
        for line in amplicon_info:
            line = line.strip()
            if not line or line[0] == '#':
                continue
            for primer_name in line.split('\t'):
                primer_amplicon_mapping[primer_name] = amplicon_id
            amplicon_id += 1

        amplicons = {}
        for primer_dat, ref_id, pool_id in cls.read_primer_bed(primer_scheme):
            name = primer_dat[2]
            try:
                mapped_id = primer_amplicon_mapping[name]
            except KeyError:
                raise ValueError(
                    f'BED file primer: "{name}" not listed in amplicon info!'
                )
            if mapped_id in amplicons:
                amplicon = amplicons[mapped_id]
                if amplicon[1] != ref_id:
                    raise ValueError('ref error')
                if amplicon[2] != pool_id:
                    raise ValueError('pool_error')
                amplicon[0].append(primer_dat)
            else:
                amplicons[mapped_id] = ([primer_dat], ref_id, pool_id)

        return cls(amplicons, scheme_name)

    @classmethod
    def read_primer_bed(cls, primer_bed):
        for record in primer_bed:
            if not record.strip() or record[0] == '#':
                continue
            fields = record.strip('\n').split('\t')

            primer_dat = (
                int(fields[1]), # start
                int(fields[2]), # end
                fields[3],      # name
                fields[5]       # strand
            )
            ref_id = fields[0]
            pool_id = fields[4]

            yield primer_dat, ref_id, pool_id

    def __init__(self, amplicons, name=None):
        self.amplicons = amplicons
        self.name = name

    def write_sanitized_bed(self, o):
        records = sorted(
            ((primer, ref, pool)
            for primers, ref, pool in self.amplicons.values()
            for primer in primers),
            # sort by ref, start, end
            key=lambda x: (x[1], x[0][0], x[0][1])
        )
        for primer_dat, ref_id, pool_id in records:
            o.write(
                f'{ref_id}\t{primer_dat[0]}\t{primer_dat[1]}\t{primer_dat[2]}\t60\t{primer_dat[3]}\n'
            )

    def write_amplicon_info(self, o, mode='full'):
        for amplicon in self.amplicons.values():
            if mode == 'full':
                names = [primer_dat[2] for primer_dat in amplicon[0]]
                o.write('\t'.join(names) + '\n')
            elif mode == 'outer':
                first = min(
                    (primer_dat for primer_dat in amplicon[0] if primer_dat[3] == '+'),
                    key=lambda x: x[0]
                )
                last = max(
                    (primer_dat for primer_dat in amplicon[0] if primer_dat[3] == '-'),
                    key=lambda x: x[1]
                )
                o.write(f'{first[2]}\t{last[2]}\n')
            elif mode == 'inner':
                first = max(
                    (primer_dat for primer_dat in amplicon[0] if primer_dat[3] == '+'),
                    key=lambda x: x[0]
                )
                last = min(
                    (primer_dat for primer_dat in amplicon[0] if primer_dat[3] == '-'),
                    key=lambda x: x[1]
                )
                o.write(f'{first[2]}\t{last[2]}\n')

    def write_insert_bed(self, o):
        for amplicon_id, amplicon in sorted(self.amplicons.items()):
            first = max(
                (primer_dat for primer_dat in amplicon[0] if primer_dat[3] == '+'),
                key=lambda x: x[0]
            )
            last = min(
                (primer_dat for primer_dat in amplicon[0] if primer_dat[3] == '-'),
                key=lambda x: x[1]
            )
            insert_start = first[1]
            insert_end = last[0]
            if self.name:
                insert_name = f'{self.name}_INSERT_{amplicon_id}'
            else:
                insert_name = f'INSERT_{amplicon_id}'
            o.write(
                f'{amplicon[1]}\t{insert_start}\t{insert_end}\t{insert_name}\t{amplicon[2]}\t+\n'
                )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Write an amplicon insert bed file for Cojac '
                    'from a BED file describing primer positions '
                    'and an (optional) amplicon info file for iVar'
    )
    parser.add_argument(
        'primer_file', type=argparse.FileType(), help='Primer file'
    )
    parser.add_argument(
        'output_file', type=argparse.FileType('w'),
        help='Output file: amplicon info file for Cojac'
    )
    parser.add_argument(
        '--ampl_file', type=argparse.FileType(), default=argparse.SUPPRESS,
        help='Amplicon info file'
    )
    args = parser.parse_args()
    print(args)
    if 'ampl_file' in args:
        scheme = Scheme.from_primers_and_amplicons(args.primer_file, args.ampl_file)
    else:
        scheme = Scheme.infer_from_primer_scheme(args.primer_file)

    scheme.write_insert_bed(args.output_file)
