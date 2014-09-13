#!/usr/bin/env python

import sys
import os
import json
import shlex
import datetime
import subprocess

def main():

    today = datetime.date.today()
    gemini_root_dir = os.environ['GEMINI_ROOT_DIR']
    params = json.loads( open( sys.argv[1] ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path']
    os.mkdir( target_directory )
    gemini_exec = os.path.join( gemini_root_dir, 'gemini', 'gemini', 'install-data.py' )
    cmd = gemini_exec + " %s %s" % (' '.join( [params['param_dict']['gerp_bp'], params['param_dict']['cadd']] ), target_directory)
    #cmd = gemini_exec + " --help > %s/foo.txt" % target_directory
    ret = subprocess.check_call( cmd, shell=True )
    data_manager_dict = { 
                'data_tables': 
                    {'gemini_databases': [ 
                            {'value': today.isoformat(), 'dbkey': 'hg19', 'name': 'GEMINI annotations (%s)' % today.isoformat(), 'path': './%s' % today.isoformat() } 
                                        ] 
                    }
                }

    #save info to json file
    with open( sys.argv[1], 'wb' ) as out:
        out.write( json.dumps( data_manager_dict ) )

if __name__ == "__main__":
    main()

