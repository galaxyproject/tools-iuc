#!/usr/bin/env python

import datetime
import json
import os
import subprocess
import sys


def main():
    today = datetime.date.today()
    params = json.loads( open( sys.argv[1] ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path']
    os.mkdir( target_directory )
    cmd = "gemini --annotation-dir %s update --dataonly %s %s" % (target_directory, params['param_dict']['gerp_bp'], params['param_dict']['cadd'] )
    subprocess.check_call( cmd, shell=True )
    data_manager_dict = {
        'data_tables': {
            'gemini_databases': [
                {'value': today.isoformat(), 'dbkey': 'hg19', 'name': 'GEMINI annotations (%s)' % today.isoformat(), 'path': './%s' % today.isoformat() }
            ]
        }
    }

    # save info to json file
    with open( sys.argv[1], 'wb' ) as out:
        out.write( json.dumps( data_manager_dict ) )


if __name__ == "__main__":
    main()
