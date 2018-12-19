#!/usr/bin/env python

import datetime
import json
import os
import subprocess
import sys

import yaml


def write_gemini_config(config, config_file):
    with open(config_file, 'w') as fo:
        yaml.dump(config, fo, allow_unicode=False, default_flow_style=False)


def main():
    today = datetime.date.today()
    params = json.loads( open( sys.argv[1] ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path']
    os.mkdir( target_directory )

    # Generate a minimal configuration file for GEMINI update
    # to instruct the tool to download the annotation data into a
    # subfolder of the target directory.
    config_file = os.path.join(target_directory, 'gemini-config.yaml')
    anno_dir = os.path.join(target_directory, 'gemini/data')
    gemini_bootstrap_config = {'annotation_dir': anno_dir}
    write_gemini_config(gemini_bootstrap_config, config_file)

    # Now gemini update can be called to download the data.
    # The GEMINI_CONFIG environment variable lets the tool discover
    # the configuration file we prepared for it.
    # Note that the tool will rewrite the file turning it into a
    # complete gemini configuration file.
    gemini_env = os.environ.copy()
    gemini_env['GEMINI_CONFIG'] = target_directory
    cmd = "gemini update --dataonly %s %s" % (
        params['param_dict']['gerp_bp'],
        params['param_dict']['cadd']
    )
    subprocess.check_call( cmd, shell=True, env=gemini_env )

    # GEMINI tool wrappers that need access to the annotation files
    # are supposed to symlink them into a gemini/data subfolder of
    # the job working directory. To have GEMINI discover them there,
    # we need to set this location as the 'annotation_dir' in the
    # configuration file.
    with open(config_file) as fi:
        config = yaml.load(fi)
    config['annotation_dir'] = 'gemini/data'
    write_gemini_config(config, config_file)

    # Finally, we prepare the metadata for the new data table record ...
    data_manager_dict = {
        'data_tables': {
            'gemini_versioned_databases': [
                {
                    'value': today.isoformat(),
                    'dbkey': 'hg19',
                    'version': params['param_dict']['gemini_db_version'],
                    'name':
                        'GEMINI annotations (%s snapshot)' % today.isoformat(),
                    'path': './%s' % today.isoformat()
                }
            ]
        }
    }

    # ... and save it to the json results file
    with open( sys.argv[1], 'wb' ) as out:
        out.write( json.dumps( data_manager_dict ) )


if __name__ == "__main__":
    main()
