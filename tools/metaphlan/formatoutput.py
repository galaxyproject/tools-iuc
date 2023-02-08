#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re
from pathlib import Path

taxo_level = {
    'k': 'kingdom',
    'p': 'phylum',
    'c': 'class',
    'o': 'order',
    'f': 'family',
    'g': 'genus',
    's': 'species',
    't': 'strains'}


def split_levels(metaphlan_output_fp, out_dp, legacy_output):
    '''
    Split default MetaPhlAn into a report for each taxonomic level

    :param metaphlan_output_fp: Path default MetaPhlAn output
    :param out_dp: Path to output directory
    :param legacy_output: Boolean for legacy output
    '''
    # prepare output files
    abund_f = {
        'k': open(out_dp / Path('kingdom'), 'w'),
        'p': open(out_dp / Path('phylum'), 'w'),
        'c': open(out_dp / Path('class'), 'w'),
        'o': open(out_dp / Path('order'), 'w'),
        'f': open(out_dp / Path('family'), 'w'),
        'g': open(out_dp / Path('genus'), 'w'),
        's': open(out_dp / Path('species'), 'w'),
        't': open(out_dp / Path('strains'), 'w')
    }
    for level in abund_f:
        abund_f[level].write("%s\t" % taxo_level[level])
        if not legacy_output:
            abund_f[level].write("%s_id\t" % taxo_level[level])
        abund_f[level].write("abundance\n")

    levels_number = len(taxo_level)

    with open(metaphlan_output_fp, 'r') as metaphlan_output_f:
        with open(out_dp / Path('all'), 'w') as all_level_f:
            # write header in all leve file
            for level in ['k', 'p', 'c', 'o', 'f', 'g', 's', 't']:
                all_level_f.write("%s\t" % taxo_level[level])
                if not legacy_output:
                    all_level_f.write("%s_id\t" % taxo_level[level])
            all_level_f.write("abundance\n")

            # parse metaphlan file
            for line in metaphlan_output_f.readlines():
                # skip headers
                if line.startswith("#"):
                    continue

                # skip UNKNOWN (v3) or UNCLASSIFIED (v4) lines in predicted taxon relative abundances
                if "UNKNOWN" in line or 'UNCLASSIFIED' in line:
                    continue

                # spit lines
                split_line = line[:-1].split('\t')
                taxo_n = split_line[0].split('|')
                if legacy_output:
                    abundance = split_line[1]
                else:
                    taxo_id = split_line[1].split('|')
                    abundance = split_line[2]

                # get taxon name and ids
                for i in range(len(taxo_n)):
                    taxo = taxo_n[i].split('__')[1]
                    taxo = taxo.replace("_", " ")
                    all_level_f.write("%s\t" % taxo)
                    if not legacy_output:
                        all_level_f.write("%s\t" % taxo_id[i])

                # if not all taxon levels
                for i in range(len(taxo_n), levels_number):
                    all_level_f.write('\t')

                all_level_f.write("%s\n" % abundance)

                # write
                last_taxo_level = taxo_n[-1].split('__')
                taxo = last_taxo_level[1].replace("_", " ")
                level = last_taxo_level[0]
                abund_f[level].write("%s\t" % taxo)
                if not legacy_output:
                    abund_f[level].write("%s\t" % taxo_id[-1])
                abund_f[level].write("%s\n" % abundance)

    # close files
    for taxo_level_f in abund_f:
        abund_f[taxo_level_f].close()


def format_for_krona(metaphlan_output_fp, krona_out_fp):
    '''
    Split default MetaPhlAn into a report for each taxonomic levKRONAel

    :param metaphlan_output_fp: Path default MetaPhlAn output
    :param krona_out: Path to output file for Krona
    '''
    re_replace = re.compile(r"\w__")
    re_bar = re.compile(r"\|")
    re_underscore = re.compile(r"_")

    with open(metaphlan_output_fp, 'r') as metaphlan_output_f:
        with open(krona_out_fp, 'w') as krona_out_f:
            for line in metaphlan_output_f.readlines():
                if "s__" in line:
                    x = line.rstrip().split('\t')
                    lineage = re.sub(re_bar, '', x[0])
                    lineage = re.sub(re_replace, '\t', lineage)
                    lineage = re.sub(re_underscore, ' ', lineage)
                    krona_out_f.write("%s\t%s\n" % (x[-1], lineage))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Format MetaPhlAn output')
    subparsers = parser.add_subparsers(dest='function')
    # split_levels
    split_levels_parser = subparsers.add_parser('split_levels', help='Split default MetaPhlAn into a report for each taxonomic level')
    split_levels_parser.add_argument('--metaphlan_output', help="Path to default MetaPhlAn output")
    split_levels_parser.add_argument('--outdir', help="Path to output directory")
    split_levels_parser.add_argument('--legacy-output', dest='legacy_output', action='store_true', help="Old MetaPhlAn2 two columns output")
    split_levels_parser.set_defaults(legacy_output=False)
    # format_for_krona
    format_for_krona_parser = subparsers.add_parser('format_for_krona', help='Split default MetaPhlAn into a report for each taxonomic level')
    format_for_krona_parser.add_argument('--metaphlan_output', help="Path to default MetaPhlAn output")
    format_for_krona_parser.add_argument('--krona_output', help="Path to Krona output directory")

    args = parser.parse_args()

    if args.function == 'split_levels':
        split_levels(
            Path(args.metaphlan_output),
            Path(args.outdir),
            args.legacy_output)
    elif args.function == 'format_for_krona':
        format_for_krona(
            Path(args.metaphlan_output),
            Path(args.krona_output))
