#!/usr/bin/env python

import csv

def process_files(primer_file, ampl_file, output_file):
    # Read the primer file
    primer_data = {}
    with open(primer_file, 'r') as primer:
        reader = csv.reader(primer, delimiter='\t')
        for row in reader:
            primer_data[row[3]] = {
                'chrom': row[0],
                'start': row[1],
                'end': row[2],
                'name': row[3],
                'pool': row[4],
                'strand': row[5]
            }

    # Process the ampl_info file and generate the output file
    with open(ampl_file, 'r') as ampl, open(output_file, 'w', newline='') as output:
        writer = csv.writer(output, delimiter='\t')

        for row in csv.reader(ampl, delimiter='\t'):
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


process_files('primers.bed', 'ampl_info.tsv', 'bed_cojac.bed')