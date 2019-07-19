#!/usr/bin/env python
# Based heavily on the HISAT2 data manager wrapper

from __future__ import print_function

import argparse
import os
import subprocess
import sys
from json import dumps, loads

DEFAULT_DATA_TABLE_NAME = "bwameth_indexes"


def build_bwameth_index(data_manager_dict, params, args):
    data_table_name = args.data_table_name
    target_directory = params[ 'output_data' ][0]['extra_files_path']
    if not os.path.exists( target_directory ):
        os.mkdir( target_directory )
    fasta_base_name = os.path.basename(args.fasta_filename)
    sym_linked_fasta_filename = os.path.join(target_directory, fasta_base_name)
    os.symlink(os.path.abspath(args.fasta_filename), sym_linked_fasta_filename)
    cmd = ['bwameth.py', 'index', sym_linked_fasta_filename]
    proc = subprocess.Popen(args=cmd, shell=False, cwd=target_directory)
    return_code = proc.wait()
    if return_code:
        print("Error building index.", file=sys.stderr)
        sys.exit( return_code )
    data_table_entry = dict(value=args.dbkey, dbkey=args.dbkey, name=args.name, path=sym_linked_fasta_filename)
    _add_data_table_entry(data_manager_dict, data_table_name, data_table_entry)


def _add_data_table_entry( data_manager_dict, data_table_name, data_table_entry ):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables'][ data_table_name ] = data_manager_dict['data_tables'].get( data_table_name, [] )
    data_manager_dict['data_tables'][ data_table_name ].append( data_table_entry )
    return data_manager_dict


def main():
    # Parse Command Line
    parser = argparse.ArgumentParser()
    parser.add_argument( '--output', default=None )
    parser.add_argument( '--fasta_filename', default=None )
    parser.add_argument( '--dbkey', default=None )
    parser.add_argument( '--name', default=None )
    parser.add_argument( '--description', default=None )
    parser.add_argument( '--data_table_name', default=DEFAULT_DATA_TABLE_NAME )
    args = parser.parse_args()

    filename = args.output
    params = loads(open(filename).read())
    data_manager_dict = {}

    if args.dbkey in [ None, '', '?' ]:
        raise Exception('"%s" is not a valid dbkey. You must specify a valid dbkey.' % (args.dbkey))

    # build the index
    build_bwameth_index(data_manager_dict, params, args)

    # save info to json file
    open(filename, 'w').write(dumps(data_manager_dict))


if __name__ == "__main__":
    main()
