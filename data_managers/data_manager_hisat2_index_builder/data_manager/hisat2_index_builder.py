#!/usr/bin/env python
# Based heavily on the Bowtie 2 data manager wrapper script by Dan Blankenberg
from __future__ import print_function

import argparse
import os
import shlex
import subprocess
import sys
from json import dumps, loads

DEFAULT_DATA_TABLE_NAME = "hisat2_indexes"


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


def build_hisat_index( data_manager_dict, options, params, sequence_id, sequence_name ):
    data_table_name = options.data_table_name or DEFAULT_DATA_TABLE_NAME
    target_directory = params[ 'output_data' ][0]['extra_files_path']
    if not os.path.exists( target_directory ):
        os.mkdir( target_directory )
    fasta_base_name = os.path.split( options.fasta_filename )[-1]
    sym_linked_fasta_filename = os.path.join( target_directory, fasta_base_name )
    os.symlink( options.fasta_filename, sym_linked_fasta_filename )
    args = [ 'hisat2-build' ]
    args.extend( shlex.split( options.indexer_options ) )
    args.extend( [ sym_linked_fasta_filename, sequence_id ] )
    threads = os.environ.get('GALAXY_SLOTS')
    if threads:
        args.extend(['-p', threads])
    proc = subprocess.Popen( args=args, shell=False, cwd=target_directory )
    return_code = proc.wait()
    if return_code:
        print("Error building index.", file=sys.stderr)
        sys.exit( return_code )
    data_table_entry = dict( value=sequence_id, dbkey=options.fasta_dbkey, name=sequence_name, path=sequence_id )
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
    parser.add_argument( '--data_table_name', dest='data_table_name', action='store', type=str, default='hisat2_indexes' )
    parser.add_argument( '--indexer_options', dest='indexer_options', action='store', type=str, default='' )
    options = parser.parse_args()

    filename = options.output

    params = loads( open( filename ).read() )
    data_manager_dict = {}

    if options.fasta_dbkey in [ None, '', '?' ]:
        raise Exception( '"%s" is not a valid dbkey. You must specify a valid dbkey.' % ( options.fasta_dbkey ) )

    sequence_id, sequence_name = get_id_name( params, dbkey=options.fasta_dbkey, fasta_description=options.fasta_description )

    # build the index
    build_hisat_index( data_manager_dict, options, params, sequence_id, sequence_name )

    # save info to json file
    open( filename, 'w' ).write( dumps( data_manager_dict ) )


if __name__ == "__main__":
    main()
