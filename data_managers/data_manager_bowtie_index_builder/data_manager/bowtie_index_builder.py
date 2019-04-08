#!/usr/bin/env python
from __future__ import print_function

import json
import optparse
import os
import subprocess
import sys
import tempfile

CHUNK_SIZE = 2**20

DEFAULT_DATA_TABLE_NAME = "bowtie_indexes"


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


def build_bowtie_index( data_manager_dict, fasta_filename, params, target_directory, dbkey, sequence_id, sequence_name, data_table_name=DEFAULT_DATA_TABLE_NAME, color_space=False ):
    # TODO: allow multiple FASTA input files
    fasta_base_name = os.path.split( fasta_filename )[-1]
    sym_linked_fasta_filename = os.path.join( target_directory, fasta_base_name )
    os.symlink( fasta_filename, sym_linked_fasta_filename )
    args = [ 'bowtie-build' ]
    if color_space:
        args.append( '-C' )
    args.append( sym_linked_fasta_filename)
    args.append( fasta_base_name )
    args.append( sym_linked_fasta_filename )
    tmp_stderr = tempfile.NamedTemporaryFile( prefix="tmp-data-manager-bowtie-index-builder-stderr" )
    proc = subprocess.Popen( args=args, shell=False, cwd=target_directory, stderr=tmp_stderr.fileno() )
    return_code = proc.wait()
    if return_code:
        tmp_stderr.flush()
        tmp_stderr.seek(0)
        print("Error building index:", file=sys.stderr)
        while True:
            chunk = tmp_stderr.read( CHUNK_SIZE )
            if not chunk:
                break
            sys.stderr.write( chunk )
        sys.exit( return_code )
    tmp_stderr.close()
    data_table_entry = dict( value=sequence_id, dbkey=dbkey, name=sequence_name, path=fasta_base_name )
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
    parser.add_option( '-n', '--data_table_name', dest='data_table_name', action='store', type="string", default=None, help='data_table_name' )
    parser.add_option( '-c', '--color_space', dest='color_space', action='store_true', default=False, help='color_space' )
    (options, args) = parser.parse_args()

    filename = args[0]

    params = json.loads( open( filename ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path']
    os.mkdir( target_directory, mode=0o755 )
    data_manager_dict = {}

    dbkey = options.fasta_dbkey

    if dbkey in [ None, '', '?' ]:
        raise Exception( '"%s" is not a valid dbkey. You must specify a valid dbkey.' % ( dbkey ) )

    sequence_id, sequence_name = get_id_name( params, dbkey=dbkey, fasta_description=options.fasta_description )

    # build the index
    build_bowtie_index( data_manager_dict, options.fasta_filename, params, target_directory, dbkey, sequence_id, sequence_name, data_table_name=options.data_table_name or DEFAULT_DATA_TABLE_NAME, color_space=options.color_space )

    # save info to json file
    open( filename, 'wb' ).write( json.dumps( data_manager_dict ) )


if __name__ == "__main__":
    main()
