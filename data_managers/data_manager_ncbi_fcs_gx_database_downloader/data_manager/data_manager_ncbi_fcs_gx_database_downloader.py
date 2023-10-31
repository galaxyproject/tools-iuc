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
    output_file: str = ''
    locfile_path: str = ''
    output_directory: str = ''
    tool_directory: str = ''
  

def main(debug=False):
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_file', required=True)

    command_line_args = parser.parse_args()

    opts = Opts()
    opts.output_file = command_line_args.output_file

    with open(opts.output_file) as f:
        params = json.load(f)

    if debug:
        with open('/var/tmp/dump.json', 'w') as ofh:
            print(json.dumps(params, indent=4), file=ofh)

    opts.tag = params['param_dict'].get('tag', '')
    opts.description = params['param_dict'].get('description', '')

    opts.default = params['param_dict'].get('default', 'false')
    if opts.default  == 'false':
        opts.default = '0'
    else:
        opts.default = '1'

    opts.phone_home = params['param_dict'].get('phone_home', 'false') 
    if opts.phone_home == 'false':
        opts.phone_home = '0'
    else:
        opts.phone_home = '1'

    opts.manifest = params['param_dict'].get('manifest', '')
    opts.local_cache_dir = params['param_dict'].get('local_cache_dir', '')
    opts.tool_cache_dir = params['param_dict'].get('tool_cache_dir', '')
    opts.output_directory = os.path.join(params['output_data'][0]['extra_files_path'], opts.tag)
    opts.tool_directory = params['param_dict']['__tool_directory__']

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



    output_dict = {
        'data_tables': {
            'ncbi_fcs_gx_databases': [
                {
                    'value': opts.tag,
                    'name': opts.description,
                    'default': opts.default,
                    'phone_home': opts.phone_home,
                    'manifest': opts.manifest,
                    'local_cache_dir': opts.local_cache_dir,
                    'tool_cache_dir': opts.tool_cache_dir,
                }
            ]
        }
    }

    with open(opts.output_file, 'w') as ofh:
        print(json.dumps(output_dict), file=ofh)

if __name__ == '__main__':
    main(debug=True)

