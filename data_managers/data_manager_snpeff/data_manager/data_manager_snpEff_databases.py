#!/usr/bin/env python
import json
import optparse
import os
import subprocess
import sys


def fetch_databases(data_manager_dict, target_directory):
    if not os.path.exists(target_directory):
        os.makedirs(target_directory)
    databases_path = os.path.join( target_directory, 'databases.out' )
    databases_output = open(databases_path, 'w')
    args = ['snpEff', 'databases']
    return_code = subprocess.call(args=args, shell=False, stdout=databases_output.fileno())
    if return_code:
        sys.exit( return_code )
    databases_output.close()
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables']['snpeffv_databases'] = data_manager_dict['data_tables'].get( 'snpeffv_databases', [] )
    data_table_entries = []
    with open(databases_path, 'r') as fh:
        for i, line in enumerate(fh):
            fields = line.split('\t')
            if len(fields) >= 2:
                genome_version = fields[0].strip()
                if genome_version.startswith("Genome") or genome_version.startswith("-"):
                    continue
                # snpeff test genome
                if genome_version == '30c2c903' or fields[1].strip() == 'TestCase' or fields[1].strip().startswith('Test_'):
                    continue
                description = fields[1].strip() + ' : ' + genome_version
                data_table_entries.append(dict(value=genome_version, name=description))
    data_manager_dict['data_tables']['snpeffv_databases'] = data_table_entries
    return data_manager_dict


def main():
    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()

    filename = args[0]

    params = json.loads( open( filename ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path']
    os.mkdir( target_directory )
    data_manager_dict = {}

    # Create Defuse Reference Data
    data_manager_dict = fetch_databases( data_manager_dict, target_directory)

    # save info to json file
    open( filename, 'w' ).write( json.dumps( data_manager_dict ) )


if __name__ == "__main__":
    main()
