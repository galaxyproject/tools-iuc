#!/usr/bin/env python

import argparse
import uuid
import os
from datetime import datetime
import json
import subprocess
import sys

def main():
    #Parse Command Line
    parser = argparse.ArgumentParser(description='Create data manager JSON.')
    parser.add_argument('--out', dest='output', action='store', help='JSON filename')
    parser.add_argument("--test", action='store_true', help="option to test the script with an lighted database")

    args = parser.parse_args()

    # the output file of a DM is a json containing args that can be used by the DM
    # most tools mainly use these args to find the extra_files_path for the DM, which can be used
    # to store the DB data
    with open(args.output) as fh:
        params = json.load(fh)


    workdir = params['output_data'][0]['extra_files_path']
    os.mkdir(workdir)

    time = datetime.utcnow().strftime("%Y-%m-%dT%H%M%SZ")
    db_value = "db_from_{0}".format(time)
    db_path = os.path.join(workdir, db_value)

    #create DB
    if args.test: #the test just checks that the pharokka download script is available
        command_args = ["install_databases.py", "-h"]   
    else:
        command_args = ["install_databases.py", "-o", db_path]

    # command_args = ["mkdir", db_path]
    proc = subprocess.Popen(args=command_args, shell=False)
    return_code = proc.wait()
    if return_code:
        print("Error downloading Pharokka database.", file=sys.stderr)
        sys.exit(return_code)

    # Update Data Manager JSON and write to file
    data_manager_entry = {
        'data_tables': {
            'pharokka_db': {
                'value': db_value,
                'dbkey': db_value,
                'name': 'Pharokka DB downloaded at {0}'.format(datetime.now()),
                'path': db_path
            }
        }
    }

    with open(os.path.join(args.output), 'w+') as fh:
        json.dump(data_manager_entry, fh, sort_keys=True)

if __name__ == '__main__':
    main()