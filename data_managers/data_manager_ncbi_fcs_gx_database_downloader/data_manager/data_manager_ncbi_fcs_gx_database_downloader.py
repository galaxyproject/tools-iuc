#!/usr/bin/env python

import argparse
import dataclasses
import json
import os
import subprocess
import sys


@dataclasses.dataclass
class Opts:
    tag: str = ''
    manifest: str = ''
    path: str = ''
    locfile_path: str = ''
    output_directory: str = ''
    output_dict: dict = dataclasses.field(default_factory=lambda: dict())


def get_options():
    opts = Opts()

    parser = argparse.ArgumentParser()
    parser.add_argument('--output_file', required=True)
    parser.add_argument('--tag', required=True)
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--path', required=True)

    command_line_args = parser.parse_args()

    opts.output_file = command_line_args.output_file
    with open(opts.output_file) as f:
        params = json.load(f)

    opts.tag = command_line_args.tag
    opts.manifest = command_line_args.manifest
    opts.path = command_line_args.path
    opts.locfile_path = os.path.join(params['param_dict'].get('__tool_directory__', ''), '..', 'test-data', 'ncbi_fcs_gx_databases.loc')
    opts.output_directory = os.path.join(params['output_data'][0]['extra_files_path'], opts.tag)

    opts.output_dict = {
        'data_tables': {
            'ncbi_fcs_gx_databases': {
                'add': [
                ],
                'remove': [
                ],
            }
        }
    }

    return opts


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

    gx_divisions = set()

    taxa_file = os.path.join(opts.output_directory, f'{opts.tag}.taxa.tsv')
    with open(taxa_file) as f:
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

    # add unknown value
    elements.append(('Unknown / Unclassified', 'unkn:unknown'))

    divisions_file = os.path.join(opts.output_directory, 'ncbi_fcs_gx_divisions.tsv')
    with open(divisions_file, 'w') as f:
        for name, value in sorted(elements):
            print(f'{value}\t{name}', file=f)


def check_locfile(opts):
    current = {}

    try:
        with open(opts.locfile_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                line = line.rstrip('\n')
                tag, manifest, path = line.split('\t', 2)
                current[tag] = (manifest, path)
    except FileNotFoundError:
        pass

    if opts.tag in current:
        current_manifest, current_path = current[opts.tag]
        if current_manifest != opts.manifest or current_path != opts.path:
            opts.output_dict['data_tables']['ncbi_fcs_gx_databases']['remove'].append({
                'value': opts.tag,
                'manifest': current_manifest,
                'name': current_path,
            })
            opts.output_dict['data_tables']['ncbi_fcs_gx_databases']['add'].append({
                'value': opts.tag,
                'manifest': opts.manifest,
                'name': opts.path,
            })
    else:
        opts.output_dict['data_tables']['ncbi_fcs_gx_databases']['add'].append({
            'value': opts.tag,
            'manifest': opts.manifest,
            'name': opts.path,
        })


def main():
    opts = get_options()
    check_locfile(opts)

    os.makedirs(opts.output_directory, exist_ok=True)

    args = [
        'sync_files.py',
        '--mft',
        opts.manifest,
        '--dir',
        opts.output_directory,
        'get',
    ]

    try:
        subprocess.run(args, capture_output=True, check=True)
    except subprocess.CalledProcessError:
        raise

    create_divisions_file(opts)

    with open(opts.output_file, 'w') as f:
        print(json.dumps(opts.output_dict, sort_keys=True, indent=2), file=f)

    sys.exit(0)


if __name__ == '__main__':
    main()
