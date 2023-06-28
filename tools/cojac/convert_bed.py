#!/usr/bin/env python

import argparse
import csv


def process_files(primer_file, ampl_file, output_file):
    primer_data = {}
    with primer_file:
        reader = csv.reader(primer_file, delimiter='\t')
        for row in reader:
            primer_data[row[3]] = {
                'chrom': row[0],
                'start': row[1],
                'end': row[2],
                'name': row[3],
                'pool': row[4],
                'strand': row[5]
            }

    with ampl_file, open(output_file, 'w', newline='') as output:
        writer = csv.writer(output, delimiter='\t')

        for row in csv.reader(ampl_file, delimiter='\t'):
            name1 = row[0]
            name2 = row[1]

            if name1 in primer_data:
                primer1 = primer_data[name1]
                primer2 = primer_data[name2]
                chrom = primer1['chrom']
                start = primer1['end']
                end = primer2['start']
                name = primer1['name'].replace('_LEFT', '').replace('_', '_INSERT_')
                score = primer1['pool'].replace('pool_', '')
                strand = '+'

                writer.writerow([chrom, start, end, name, score, strand])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Write an amplicon info file for Cojac '
                    'from a BED file describing primer positions '
                    'and amplicon info file for iVar'
    )
    parser.add_argument(
        'primer_file', type=argparse.FileType(), help='Primer file'
    )
    parser.add_argument(
        'ampl_file', type=argparse.FileType(), help='Amplicon info file'
    )
    parser.add_argument(
        'output_file', type=argparse.FileType('w'), help='Output file: amplicon info file for Cojac'
    )
    args = parser.parse_args()

    process_files(args.primer_file, args.ampl_file, args.output_file.name)
