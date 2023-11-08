#!/usr/bin/env python

import argparse
import dataclasses
import json
import os
import subprocess


@dataclasses.dataclass
class Opts:
    tag: str = ''
    description: str = ''
    default: str = ''
    phone_home: str = ''
    manifest: str = ''
    local_cache_dir: str = ''
    tool_cache_dir: str = ''
    locfile_path: str = ''
    output_directory: str = ''
    output_dict: dict = dataclasses.field(default_factory=lambda: dict())


def get_options():
    opts = Opts()

    parser = argparse.ArgumentParser()
    parser.add_argument('--output_file', required=True)
    command_line_args = parser.parse_args()

    opts.output_file = command_line_args.output_file
    with open(opts.output_file) as f:
        params = json.load(f)

    opts.tag = params['param_dict'].get('tag', '')
    opts.description = params['param_dict'].get('description', '')
    opts.default = '1' if params['param_dict'].get('default', 'false') == 'true' else '0'
    opts.phone_home = '1' if params['param_dict'].get('phone_home', 'false') == 'true' else '0'
    opts.manifest = params['param_dict'].get('manifest', '')
    opts.local_cache_dir = params['param_dict'].get('local_cache_dir', '')
    opts.tool_cache_dir = params['param_dict'].get('tool_cache_dir', '')
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

    divisions_file = os.path.join(opts.output_directory, 'ncbi_fcs_gx_divisions.tsv')
    with open(divisions_file, 'w') as f:
        print('## NCBI FCS GX Divisions', file=f)
        print('#', file=f)
        print('#value\tname', file=f)
        for name, value in sorted(elements):
            print(f'{value}\t{name}', file=f)


def check_locfile(opts):
    current = {}
    have_default = False
    default_tag = ''

    try:
        with open(opts.locfile_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                line = line.rstrip('\n')
                tag, description, default, phone_home, manifest, local_cache_dir, tool_cache_dir = line.split('\t', 6)
                if default == '1' and not have_default:
                    have_default = True
                    default_tag = tag
                current[tag] = (description, default, phone_home, manifest, local_cache_dir, tool_cache_dir)
    except FileNotFoundError:
        pass

    if opts.tag in current:
        description, default, phone_home, manifest, local_cache_dir, tool_cache_dir = current[opts.tag]

        if have_default and default == '1':
            have_default = False

        opts.output_dict['data_tables']['ncbi_fcs_gx_databases']['remove'].append({
            'value': opts.tag,
            'name': description,
            'default': default,
            'phone_home': phone_home,
            'manifest': manifest,
            'local_cache_dir': local_cache_dir,
            'tool_cache_dir': tool_cache_dir,
        })

    if opts.default == '0':
        if have_default:
            description, default, phone_home, manifest, local_cache_dir, tool_cache_dir = current[default_tag]
            for operation in ['remove', 'add']:
                opts.output_dict['data_tables']['ncbi_fcs_gx_databases'][operation].append({
                    'value': default_tag,
                    'name': description,
                    'default': default,
                    'phone_home': phone_home,
                    'manifest': manifest,
                    'local_cache_dir': local_cache_dir,
                    'tool_cache_dir': tool_cache_dir,
                })
        else:
            opts.default = '1'

    opts.output_dict['data_tables']['ncbi_fcs_gx_databases']['add'].append({
        'value': opts.tag,
        'name': opts.description,
        'default': opts.default,
        'phone_home': opts.phone_home,
        'manifest': opts.manifest,
        'local_cache_dir': opts.local_cache_dir,
        'tool_cache_dir': opts.tool_cache_dir,
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


if __name__ == '__main__':
    main()
