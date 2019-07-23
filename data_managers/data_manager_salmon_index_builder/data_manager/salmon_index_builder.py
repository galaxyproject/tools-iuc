#!/usr/bin/env python
# Based heavily on the kallisto data manager wrapper script by iuc
from __future__ import print_function

import argparse
import os
import subprocess
import sys
from json import dumps, loads

DEFAULT_DATA_TABLE_NAME = "salmon_indexes_versioned"


def get_id_name( params, dbkey, fasta_description=None):
    # TODO: ensure sequence_id is unique and does not already appear in location file
    sequence_id = params['param_dict']['sequence_id']
    if not sequence_id:
        sequence_id = dbkey

    sequence_name = params['param_dict']['sequence_name']
    if not sequence_name:
        sequence_name = fasta_description
        if not sequence_name:
            sequence_name = dbkey
    return sequence_id, sequence_name


def build_salmon_index( data_manager_dict, options, params, sequence_id, sequence_name ):
    data_table_name = options.data_table_name or DEFAULT_DATA_TABLE_NAME
    target_directory = params[ 'output_data' ][0]['extra_files_path']
    if not os.path.exists( target_directory ):
        os.mkdir( target_directory )
    args = [ 'salmon', 'index' ]
    if options.kmer_size != '':
        args.append('-k')
        args.append(options.kmer_size)
    args.extend( [ '-t', options.fasta_filename, '-i', target_directory ] )
    proc = subprocess.Popen( args=args, shell=False)
    return_code = proc.wait()
    if return_code:
        print("Error building index.", file=sys.stderr)
        sys.exit( return_code )
    data_table_entry = dict( value=sequence_id, dbkey=options.fasta_dbkey, name=sequence_name, path=sequence_id, version=options.index_version )
    _add_data_table_entry( data_manager_dict, data_table_name, data_table_entry )


def _add_data_table_entry( data_manager_dict, data_table_name, data_table_entry ):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables'][ data_table_name ] = data_manager_dict['data_tables'].get( data_table_name, [] )
    data_manager_dict['data_tables'][ data_table_name ].append( data_table_entry )
    return data_manager_dict


def main():
    # Parse Command Line
    parser = argparse.ArgumentParser()
    parser.add_argument( '--output', dest='output', action='store', type=str, default=None )
    parser.add_argument( '--fasta_filename', dest='fasta_filename', action='store', type=str, default=None )
    parser.add_argument( '--fasta_dbkey', dest='fasta_dbkey', action='store', type=str, default=None )
    parser.add_argument( '--fasta_description', dest='fasta_description', action='store', type=str, default=None )
    parser.add_argument( '--data_table_name', dest='data_table_name', action='store', type=str, default='salmon_indexes' )
    parser.add_argument( '-v', '--index_version', dest='index_version', action='store', type=str, help='Use IndexVersion attribute from header.json' )
    parser.add_argument( '-k', '--kmer_size', dest='kmer_size', action='store', type=str, help='kmer_size' )
    options = parser.parse_args()

    filename = options.output

    params = loads( open( filename ).read() )
    data_manager_dict = {}

    if options.fasta_dbkey in [ None, '', '?' ]:
        raise Exception( '"%s" is not a valid dbkey. You must specify a valid dbkey.' % ( options.fasta_dbkey ) )

    sequence_id, sequence_name = get_id_name( params, dbkey=options.fasta_dbkey, fasta_description=options.fasta_description )
    # build the index
    build_salmon_index( data_manager_dict, options, params, sequence_id, sequence_name )

    # save info to json file
    open( filename, 'w' ).write( dumps( data_manager_dict ) )


if __name__ == "__main__":
    main()
