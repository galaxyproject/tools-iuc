#!/usr/bin/env python

import datetime
import json
import os
import subprocess
import sys

import yaml


def main():
    today = datetime.date.today()
    params = json.loads( open( sys.argv[1] ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path']
    os.mkdir( target_directory )
    # The target_directory needs to be specified twice for the following
    # invocation of gemini.
    # In essence, the GEMINI_CONFIG environment variable makes gemini store
    # its yaml configuration file in that directory, while the
    # --annotation-dir argument makes it write the same path into the yaml
    # file, which is then used for determining where the actual annotation
    # files should be stored.
    gemini_env = os.environ.copy()
    gemini_env['GEMINI_CONFIG'] = target_directory
    cmd = "gemini --annotation-dir %s update --dataonly %s %s" % (
        target_directory,
        params['param_dict']['gerp_bp'],
        params['param_dict']['cadd']
    )
    subprocess.check_call( cmd, shell=True, env=gemini_env )

    # modify the newly created gemini config file to contain a relative
    # annotation dir path, which will be interpreted as relative to
    # the job working directory at runtime by any gemini tool
    config_file = os.path.join(target_directory, 'gemini-config.yaml')
    with open(config_file) as fi:
        config = yaml.load(fi)
    config['annotation_dir'] = 'gemini/data'
    with open(config_file, 'w') as fo:
        yaml.dump(config, fo, allow_unicode=False, default_flow_style=False)

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

    # save info to json file
    with open( sys.argv[1], 'wb' ) as out:
        out.write( json.dumps( data_manager_dict ) )


if __name__ == "__main__":
    main()
