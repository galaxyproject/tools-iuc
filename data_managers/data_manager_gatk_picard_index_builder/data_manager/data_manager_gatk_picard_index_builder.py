#!/usr/bin/env python
# Dave B.
# Uses fasta sorting functions written by Dan Blankenberg.

import json
import optparse
import os
import shutil
import subprocess
import sys
import tempfile

CHUNK_SIZE = 2**20
DEFAULT_DATA_TABLE_NAME = "fasta_indexes"


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


def build_picard_index( data_manager_dict, fasta_filename, target_directory, dbkey, sequence_id, sequence_name, jar, data_table_name=DEFAULT_DATA_TABLE_NAME ):
    fasta_base_name = os.path.split( fasta_filename )[-1]
    gatk_sorted_fasta_filename = os.path.join( target_directory, fasta_base_name )
    shutil.copy( fasta_filename, gatk_sorted_fasta_filename )
    _sort_fasta_gatk( gatk_sorted_fasta_filename )
    sam_index_filename = '%s.fai' % gatk_sorted_fasta_filename
    if not os.path.exists( sam_index_filename ):
        sam_command = [ 'samtools', 'faidx', gatk_sorted_fasta_filename ]
        _run_command( sam_command, target_directory )
    args = [ 'java', '-jar', jar, 'R=%s' % gatk_sorted_fasta_filename, 'O=%s.dict' % sequence_id ]
    _run_command( args, target_directory )
    data_table_entry = dict( value=sequence_id, dbkey=dbkey, name=sequence_name, path=fasta_base_name )
    _add_data_table_entry( data_manager_dict, data_table_name, data_table_entry )


def _run_command( command, target_directory ):
    tmp_stderr = tempfile.NamedTemporaryFile( prefix="tmp-data-manager-gatk_picard_index_builder-stderr" )
    proc = subprocess.Popen( args=command, shell=False, cwd=target_directory, stderr=tmp_stderr.fileno() )
    return_code = proc.wait()
    if return_code:
        tmp_stderr.flush()
        tmp_stderr.seek( 0 )
        sys.stderr.write( "Error building index:\n" )
        while True:
            chunk = tmp_stderr.read( CHUNK_SIZE )
            if not chunk:
                break
            sys.stderr.write( chunk )
        sys.exit( return_code )
    tmp_stderr.close()


def _add_data_table_entry( data_manager_dict, data_table_name, data_table_entry ):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables'][ data_table_name ] = data_manager_dict['data_tables'].get( data_table_name, [] )
    data_manager_dict['data_tables'][ data_table_name ].append( data_table_entry )
    return data_manager_dict


def _move_and_index_fasta_for_sorting( fasta_filename ):
    unsorted_filename = tempfile.NamedTemporaryFile().name
    shutil.move( fasta_filename, unsorted_filename )
    fasta_offsets = {}
    unsorted_fh = open( unsorted_filename )
    while True:
        offset = unsorted_fh.tell()
        line = unsorted_fh.readline()
        if not line:
            break
        if line.startswith( ">" ):
            line = line.split( None, 1 )[0][1:]
            fasta_offsets[ line ] = offset
    unsorted_fh.close()
    current_order = [x[1] for x in sorted( ( x[1], x[0] ) for x in fasta_offsets.items() )]
    return ( unsorted_filename, fasta_offsets, current_order )


def _write_sorted_fasta( sorted_names, fasta_offsets, sorted_fasta_filename, unsorted_fasta_filename ):
    unsorted_fh = open( unsorted_fasta_filename )
    sorted_fh = open( sorted_fasta_filename, 'wb+' )

    for name in sorted_names:
        offset = fasta_offsets[ name ]
        unsorted_fh.seek( offset )
        sorted_fh.write( unsorted_fh.readline() )
        while True:
            line = unsorted_fh.readline()
            if not line or line.startswith( ">" ):
                break
            sorted_fh.write( line )
    unsorted_fh.close()
    sorted_fh.close()


def _int_to_roman( integer ):
    if not isinstance( integer, int ):
        raise TypeError("expected integer, got %s" % type( integer ))
    if not 0 < integer < 4000:
        raise ValueError("Argument must be between 1 and 3999, got %s" % str( integer ))
    ints = ( 1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1 )
    nums = ( 'M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I' )
    result = ""
    for i in range( len( ints ) ):
        count = int( integer / ints[ i ] )
        result += nums[ i ] * count
        integer -= ints[ i ] * count
    return result


def _sort_fasta_gatk( fasta_filename ):
    ( unsorted_filename, fasta_offsets, current_order ) = _move_and_index_fasta_for_sorting( fasta_filename )
    sorted_names = list(map( str, range( 1, 100 ) )) + list(map( _int_to_roman, range( 1, 100 ) )) + [ 'X', 'Y', 'M' ]
    # detect if we have chrN, or just N
    has_chr = False
    for chrom in sorted_names:
        if "chr%s" % chrom in current_order:
            has_chr = True
            break

    if has_chr:
        sorted_names = ["chr%s" % x for x in sorted_names]
    else:
        sorted_names.insert( 0, "MT" )
    sorted_names.extend( [ "%s_random" % x for x in sorted_names ] )

    existing_sorted_names = []
    for name in sorted_names:
        # Append each chromosome only once.
        if name in current_order and name not in existing_sorted_names:
            existing_sorted_names.append( name )
    for name in current_order:
        # TODO: confirm that non-canonical names do not need to be sorted specially
        if name not in existing_sorted_names:
            existing_sorted_names.append( name )

    if existing_sorted_names == current_order:
        shutil.move( unsorted_filename, fasta_filename )
    else:
        _write_sorted_fasta( existing_sorted_names, fasta_offsets, fasta_filename, unsorted_filename )


def main():
    parser = optparse.OptionParser()
    parser.add_option( '-f', '--fasta_filename', dest='fasta_filename', action='store', type="string", default=None, help='fasta_filename' )
    parser.add_option( '-d', '--fasta_dbkey', dest='fasta_dbkey', action='store', type="string", default=None, help='fasta_dbkey' )
    parser.add_option( '-t', '--fasta_description', dest='fasta_description', action='store', type="string", default=None, help='fasta_description' )
    parser.add_option( '-n', '--data_table_name', dest='data_table_name', action='store', type="string", default=None, help='data_table_name' )
    parser.add_option( '-j', '--jar', dest='jar', action='store', type="string", default=None, help='GATK .jar file' )
    (options, args) = parser.parse_args()

    filename = args[0]

    params = json.loads( open( filename ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path']
    os.mkdir( target_directory )
    data_manager_dict = {}

    if options.fasta_dbkey in [ None, '', '?' ]:
        raise Exception( '"%s" is not a valid dbkey. You must specify a valid dbkey.' % ( options.fasta_dbkey ) )

    sequence_id, sequence_name = get_id_name( params, dbkey=options.fasta_dbkey, fasta_description=options.fasta_description )

    # build the index
    build_picard_index( data_manager_dict,
                        options.fasta_filename,
                        target_directory,
                        options.fasta_dbkey,
                        sequence_id,
                        sequence_name,
                        options.jar,
                        data_table_name=options.data_table_name or DEFAULT_DATA_TABLE_NAME )

    # save info to json file
    open( filename, 'wb' ).write( json.dumps( data_manager_dict ) )


if __name__ == "__main__":
    main()
