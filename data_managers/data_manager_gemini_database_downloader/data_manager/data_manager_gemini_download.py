#!/usr/bin/env python2

# IMPORTANT: This will run using Python 2 still!

import datetime
import json
import os
import subprocess
import sys

import yaml


def write_gemini_config(config, config_file):
    with open(config_file, 'w') as fo:
        yaml.dump(config, fo, allow_unicode=False, default_flow_style=False)


def load_gemini_config(config_file):
    with open(config_file) as fi:
        return yaml.load(fi)


def main():
    today = datetime.date.today()
    with open(sys.argv[1]) as fh:
        params = json.load(fh)
    target_directory = params['output_data'][0]['extra_files_path']
    os.mkdir(target_directory)

    # Prepare the metadata for the new data table record

    # The name of the database should reflect whether it was built with or
    # without the optional GERP-bp data, the CADD scores, or both.
    # This builds up the correpsonding part of the name:
    anno_extras = []
    if params['param_dict']['gerp_bp']:
        anno_extras.append('GERP')
    if params['param_dict']['cadd']:
        anno_extras.append('CADD')
    if anno_extras:
        anno_desc = ' w/ ' + ' & '.join(anno_extras)
    else:
        anno_desc = ''

    data_manager_dict = {
        'data_tables': {
            'gemini_versioned_databases': [
                {
                    'value': today.isoformat(),
                    'dbkey': 'hg19',
                    'version': params['param_dict']['gemini_db_version'],
                    'name':
                        'GEMINI annotations%s (%s snapshot)' % (
                            anno_desc, today.isoformat()
                    ),
                    'path': './%s' % today.isoformat()
                }
            ]
        }
    }

    # Save the data table metadata to the json results file
    with open(sys.argv[1], 'w') as fh:
        json.dump(data_manager_dict, fh, sort_keys=True)

    # Generate a minimal configuration file for GEMINI update
    # to instruct the tool to download the annotation data into a
    # subfolder of the target directory.
    config_file = os.path.join(target_directory, 'gemini-config.yaml')
    anno_dir = os.path.join(target_directory, 'gemini/data')
    gemini_bootstrap_config = {'annotation_dir': anno_dir}
    write_gemini_config(gemini_bootstrap_config, config_file)

    # Verify that we can read the config_file just created as we need to do so
    # after the data download has finished and it is very annoying to have this
    # fail after dozens of Gbs of data have been downloaded
    config = load_gemini_config(config_file)

    # Now gemini update can be called to download the data.
    # The GEMINI_CONFIG environment variable lets the tool discover
    # the configuration file we prepared for it.
    # Note that the tool will rewrite the file turning it into a
    # complete gemini configuration file.
    gemini_env = os.environ.copy()
    gemini_env['GEMINI_CONFIG'] = target_directory
    cmd = ['gemini', 'update', '--dataonly']
    if params['param_dict']['gerp_bp']:
        cmd += ['--extra', 'gerp_bp']
    if params['param_dict']['cadd']:
        cmd += ['--extra', 'cadd_score']

    if not params['param_dict']['test_data_manager']:
        # This is not a test => Going to embark on a massive download now
        subprocess.check_call(cmd, env=gemini_env)

    # GEMINI tool wrappers that need access to the annotation files
    # are supposed to symlink them into a gemini/data subfolder of
    # the job working directory. To have GEMINI discover them there,
    # we need to set this location as the 'annotation_dir' in the
    # configuration file.
    config = load_gemini_config(config_file)
    config['annotation_dir'] = 'gemini/data'
    write_gemini_config(config, config_file)


if __name__ == "__main__":
    main()
