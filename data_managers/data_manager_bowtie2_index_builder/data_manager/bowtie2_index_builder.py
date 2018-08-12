#!/usr/bin/env python
# Dan Blankenberg
from __future__ import print_function

import optparse
import os
import subprocess
import sys
from json import dumps, loads

DEFAULT_DATA_TABLE_NAMES = ["bowtie2_indexes"]


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


def build_bowtie2_index( data_manager_dict, fasta_filename, params, target_directory, dbkey, sequence_id, sequence_name, data_table_names=DEFAULT_DATA_TABLE_NAMES ):
    # TODO: allow multiple FASTA input files
    fasta_base_name = os.path.split( fasta_filename )[-1]
    sym_linked_fasta_filename = os.path.join( target_directory, fasta_base_name )
    os.symlink( fasta_filename, sym_linked_fasta_filename )
    args = [ 'bowtie2-build', sym_linked_fasta_filename, sequence_id ]
    proc = subprocess.Popen( args=args, shell=False, cwd=target_directory )
    return_code = proc.wait()
    if return_code:
        print("Error building index.", file=sys.stderr)
        sys.exit( return_code )
    data_table_entry = dict( value=sequence_id, dbkey=dbkey, name=sequence_name, path=sequence_id )
    for data_table_name in data_table_names:
        _add_data_table_entry( data_manager_dict, data_table_name, data_table_entry )


def _add_data_table_entry( data_manager_dict, data_table_name, data_table_entry ):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables'][ data_table_name ] = data_manager_dict['data_tables'].get( data_table_name, [] )
    data_manager_dict['data_tables'][ data_table_name ].append( data_table_entry )
    return data_manager_dict


def main():
    parser = optparse.OptionParser()
    parser.add_option( '-f', '--fasta_filename', dest='fasta_filename', action='store', type="string", default=None, help='fasta_filename' )
    parser.add_option( '-d', '--fasta_dbkey', dest='fasta_dbkey', action='store', type="string", default=None, help='fasta_dbkey' )
    parser.add_option( '-t', '--fasta_description', dest='fasta_description', action='store', type="string", default=None, help='fasta_description' )
    parser.add_option( '-n', '--data_table_name', dest='data_table_name', action='append', type="string", default=None, help='data_table_name' )
    (options, args) = parser.parse_args()

    filename = args[0]

    params = loads( open( filename ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path']
    os.mkdir( target_directory )
    data_manager_dict = {}

    dbkey = options.fasta_dbkey

    if dbkey in [ None, '', '?' ]:
        raise Exception( '"%s" is not a valid dbkey. You must specify a valid dbkey.' % ( dbkey ) )

    sequence_id, sequence_name = get_id_name( params, dbkey=dbkey, fasta_description=options.fasta_description )

    # build the index
    build_bowtie2_index( data_manager_dict, options.fasta_filename, params, target_directory, dbkey, sequence_id, sequence_name, data_table_names=options.data_table_name or DEFAULT_DATA_TABLE_NAMES )

    # save info to json file
    open( filename, 'w' ).write( dumps( data_manager_dict ) )


if __name__ == "__main__":
    main()
