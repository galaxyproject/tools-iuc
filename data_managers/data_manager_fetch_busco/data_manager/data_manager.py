#!/usr/bin/env python
#
# Data manager for reference data for the 'BUSCO' Galaxy tools
import argparse
import json
import os
import shutil
import subprocess
import datetime

def main(args):
    workdir = os.path.join(os.getcwd(), 'busco_downloads')
    cmd = "busco --download %s" % args.database
    subprocess.check_call(cmd, shell=True)
    data_manager_entry = {}
    data_manager_entry['value'] = args.name.lower()
    data_manager_entry['name'] = args.name
    data_manager_entry['version'] = args.version
    data_manager_entry['path'] = '.'
    data_manager_json = dict(data_tables=dict(busco=data_manager_entry))
    with open(args.json) as fh:
        params = json.load(fh)
    target_directory = params['output_data'][0]['extra_files_path']
    os.mkdir(target_directory)
    output_path = os.path.abspath(os.path.join(os.getcwd(), 'busco_downloads'))
    for filename in os.listdir(workdir):
        shutil.move(os.path.join(output_path, filename), target_directory)
    with open(args.json, 'w') as fh:
        json.dump(data_manager_json, fh, sort_keys=True)

if __name__ == "__main__":

    # Read command line
    parser = argparse.ArgumentParser(description='Download BUSCO database')
    parser.add_argument('--database', help="Database name")
    parser.add_argument('--name', default=str(datetime.date.today()), help='Data table entry unique ID')
    parser.add_argument('--version', help="BUSCO version")
    parser.add_argument('--json', help="Path to JSON file")
    args = parser.parse_args()
    
    main(args)

