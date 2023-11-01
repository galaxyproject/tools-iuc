#!/usr/bin/env python

import argparse
import json
import os
import sys
import subprocess
import dataclasses

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


def check_locfile(opts):
    current = {}
    have_default = False

    try:
        with open(opts.locfile_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                line = line.rstrip('\n')
                tag, description, default, phone_home, manifest, local_cache_dir, tool_cache_dir = line.split('\t', 6)
                if default == '1' and not have_default:
                    have_default = True
                current[tag] = (description, default, phone_home, manifest, local_cache_dir, tool_cache_dir)
    except:
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

    if not have_default and opts.default == '0':
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

def main(debug=False):
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
        p = subprocess.run(args, capture_output=True, check=True)
    except subprocess.CalledProcessError as e:
        print ('error')
        print (f'path={os.environ["PATH"]}')
        print (f'returncode={e.returncode}')
        print (f'cmd={e.cmd}')
        print (f'output={e.output}')
        print (f'stdout={e.stdout}')
        print (f'stderr={e.stderr}')
        raise

    with open(opts.output_file, 'w') as f:
        print(json.dumps(opts.output_dict, sort_keys=True, indent=2), file=f)


if __name__ == '__main__':
    main(debug=True)

