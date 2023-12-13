#!/usr/bin/env python

import argparse
import json
import os
import subprocess


def main():
    opts = parse_args()
    sync_files(opts)
    create_divisions_file(opts)

    output_dict = {
        'data_tables': {
            'ncbi_fcs_gx_databases': {
                'add': [
                    {
                        'value': opts.tag,
                        'source_manifest': opts.source_manifest,
                        'name': opts.output_dir
                    }
                ]
            },
            'ncbi_fcs_gx_divisions': {
                'add': [
                    {
                        'value': opts.tag,
                        'name': opts.output_dir
                    }
                ]
            }
        }
    }

    with open(opts.output_file, 'w') as f:
        print(json.dumps(output_dict, sort_keys=True, indent=2), file=f)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tag', required=True)
    parser.add_argument('--source_manifest', required=True)
    parser.add_argument('--output_file', required=True)
    parser.add_argument('--output_dir', required=True)

    return parser.parse_args()


def sync_files(opts):
    db_dir = os.path.join(opts.output_dir, 'db')
    os.makedirs(db_dir, exist_ok=True)

    args = [
        'sync_files.py',
        '--mft',
        opts.source_manifest,
        '--dir',
        db_dir,
        'get'
    ]

    try:
        subprocess.run(args, capture_output=True, check=True)
    except subprocess.CalledProcessError:
        raise


def create_divisions_file(opts):
    top_level_description = {
        'anml': 'Animals (Metazoa)',
        'arch': 'Archaea',
        'fung': 'Fungi',
        'plnt': 'Plants (Viridiplantae)',
        'prok': 'Bacteria',
        'prst': 'Protists (other Eukaryota)',
        'synt': 'Synthetic',
        'virs': 'Virus',
    }

    manifest_filename = os.path.basename(opts.source_manifest)
    assert manifest_filename.lower().endswith('.manifest'), 'source_manifest does not end with ".manifest"'
    manifest_tag = manifest_filename[:-9]

    gx_divisions = set()

    taxa_pathname = os.path.join(opts.output_dir, 'db', f'{manifest_tag}.taxa.tsv')
    with open(taxa_pathname) as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.rstrip('\n')
            tax_id, species, common_name, blast_div, div = line.split('\t', 4)
            gx_divisions.add(div)

    elements = []
    for division in gx_divisions:
        top, rest = division.split(':', 1)
        description = f'{top_level_description[top]} - {rest}'
        elements.append((description, division))

    # add element to support unknown/unclassified samples
    elements.append(('Unknown / Unclassified', 'unkn:unknown'))

    divisions_file = os.path.join(opts.output_dir, 'ncbi_fcs_gx_divisions.tsv')
    with open(divisions_file, 'w') as f:
        for name, value in sorted(elements):
            print(f'{value}\t{name}', file=f)


if __name__ == '__main__':
    main()
