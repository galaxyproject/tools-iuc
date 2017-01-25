#!/usr/bin/env python
#Dan Blankenberg

import sys
import os
import tempfile
import shutil
import optparse
import urllib2
#import uuid
from ftplib import FTP
import tarfile
import zipfile
import gzip
import bz2

from json import loads, dumps


CHUNK_SIZE = 2**20 #1mb

def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)
    
def get_dbkey_id_name( params, dbkey_description=None ):
    dbkey = params['param_dict']['dbkey']
    #TODO: ensure sequence_id is unique and does not already appear in location file
    sequence_id = params['param_dict']['sequence_id']
    if not sequence_id:
        sequence_id = dbkey #uuid.uuid4() generate and use an uuid instead?
    
    sequence_name = params['param_dict']['sequence_name']
    if not sequence_name:
        sequence_name = dbkey_description
        if not sequence_name:
            sequence_name = dbkey
    return dbkey, sequence_id, sequence_name

def _get_files_in_ftp_path( ftp, path ):
    path_contents = []
    ftp.retrlines( 'MLSD %s' % ( path ), path_contents.append )
    return [ line.split( ';' )[ -1 ].lstrip() for line in path_contents ]

def _get_stream_readers_for_tar( file_obj, tmp_dir ):
    fasta_tar = tarfile.open( fileobj=file_obj, mode='r:*' )
    return filter( lambda x: x is not None, [ fasta_tar.extractfile( member ) for member in fasta_tar.getmembers() ] )

def _get_stream_readers_for_zip( file_obj, tmp_dir ):
    fasta_zip = zipfile.ZipFile( file_obj, 'r' )
    rval = []
    for member in fasta_zip.namelist():
        fasta_zip.extract( member, tmp_dir )
        rval.append( open( os.path.join( tmp_dir, member ), 'rb' ) )
    return rval

def _get_stream_readers_for_gzip( file_obj, tmp_dir ):
    return [ gzip.GzipFile( fileobj=file_obj, mode='rb' ) ]

def _get_stream_readers_for_bz2( file_obj, tmp_dir ):
    return [ bz2.BZ2File( file_obj.name, 'rb' ) ]

def sort_fasta( fasta_filename, sort_method, params ):
    if sort_method is None:
        return
    assert sort_method in SORTING_METHODS, ValueError( "%s is not a valid sorting option." % sort_method )
    return SORTING_METHODS[ sort_method ]( fasta_filename, params )

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
    current_order = map( lambda x: x[1], sorted( map( lambda x: ( x[1], x[0] ), fasta_offsets.items() ) ) )
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

def _sort_fasta_as_is( fasta_filename, params ):
    return

def _sort_fasta_lexicographical( fasta_filename, params ):
    ( unsorted_filename, fasta_offsets, current_order ) = _move_and_index_fasta_for_sorting( fasta_filename )
    sorted_names = sorted( fasta_offsets.keys() )
    if sorted_names == current_order:
        shutil.move( unsorted_filename, fasta_filename )
    else:
        _write_sorted_fasta( sorted_names, fasta_offsets, fasta_filename, unsorted_filename )    

def _sort_fasta_gatk( fasta_filename, params ):
    #This method was added by reviewer request.
    ( unsorted_filename, fasta_offsets, current_order ) = _move_and_index_fasta_for_sorting( fasta_filename )
    sorted_names = map( str, range( 1, 23 ) ) + [ 'X', 'Y' ]
    #detect if we have chrN, or just N
    has_chr = False
    for chrom in sorted_names:
        if "chr%s" % chrom in current_order:
            has_chr = True
            break
    
    if has_chr:
        sorted_names = map( lambda x: "chr%s" % x, sorted_names)
        sorted_names.insert( 0, "chrM" )
    else:
        sorted_names.insert( 0, "MT" )
    sorted_names.extend( map( lambda x: "%s_random" % x, sorted_names ) )
    
    existing_sorted_names = []
    for name in sorted_names:
        if name in current_order:
            existing_sorted_names.append( name )
    for name in current_order:
        #TODO: confirm that non-canonical names do not need to be sorted specially
        if name not in existing_sorted_names:
            existing_sorted_names.append( name )
    
    if existing_sorted_names == current_order:
        shutil.move( unsorted_filename, fasta_filename )
    else:
        _write_sorted_fasta( existing_sorted_names, fasta_offsets, fasta_filename, unsorted_filename )

def _sort_fasta_custom( fasta_filename, params ):
    ( unsorted_filename, fasta_offsets, current_order ) = _move_and_index_fasta_for_sorting( fasta_filename )
    sorted_names = []
    for id_repeat in params['param_dict']['sorting']['sequence_identifiers']:
        sorted_names.append( id_repeat[ 'identifier' ] )
    handle_not_listed = params['param_dict']['sorting']['handle_not_listed']['handle_not_listed_selector']
    if handle_not_listed.startswith( 'keep' ):
        add_list = []
        for name in current_order:
            if name not in sorted_names:
                add_list.append( name )
        if add_list:
            if handle_not_listed == 'keep_append':
                sorted_names.extend( add_list )
            else:
                add_list.extend( sorted_names )
                sorted_names = add_list
    if sorted_names == current_order:
        shutil.move( unsorted_filename, fasta_filename )
    else:
        _write_sorted_fasta( sorted_names, fasta_offsets, fasta_filename, unsorted_filename )

def download_from_ucsc( data_manager_dict, params, target_directory, dbkey, sequence_id, sequence_name ):
    UCSC_FTP_SERVER = 'hgdownload.cse.ucsc.edu'
    UCSC_DOWNLOAD_PATH = '/goldenPath/%s/bigZips/'
    COMPRESSED_EXTENSIONS = [ ( '.tar.gz', _get_stream_readers_for_tar ), ( '.tar.bz2', _get_stream_readers_for_tar ), ( '.zip', _get_stream_readers_for_zip ), ( '.fa.gz', _get_stream_readers_for_gzip ), ( '.fa.bz2', _get_stream_readers_for_bz2 ) ]
    
    email = params['param_dict']['__user_email__']
    if not email:
        email = 'anonymous@example.com'

    ucsc_dbkey = params['param_dict']['reference_source']['requested_dbkey'] or dbkey
    UCSC_CHROM_FA_FILENAMES = [ '%s.chromFa' % ucsc_dbkey, 'chromFa', ucsc_dbkey ]
    
    ftp = FTP( UCSC_FTP_SERVER )
    ftp.login( 'anonymous', email )
    
    ucsc_path = UCSC_DOWNLOAD_PATH % ucsc_dbkey
    path_contents = _get_files_in_ftp_path( ftp, ucsc_path )
    
    ucsc_file_name = None
    get_stream_reader = None
    ext = None
    ucsc_chrom_fa_filename = None
    for ucsc_chrom_fa_filename in UCSC_CHROM_FA_FILENAMES:
        for ext, get_stream_reader in COMPRESSED_EXTENSIONS:
            if "%s%s" % ( ucsc_chrom_fa_filename, ext ) in path_contents:
                ucsc_file_name = "%s%s%s" % ( ucsc_path, ucsc_chrom_fa_filename, ext )
                break
        if ucsc_file_name:
            break
    
    if not ucsc_file_name:
        raise Exception( 'Unable to determine filename for UCSC Genome for %s: %s' % ( ucsc_dbkey, path_contents ) )
    
    
    tmp_dir = tempfile.mkdtemp( prefix='tmp-data-manager-ucsc-' )
    ucsc_fasta_filename = os.path.join( tmp_dir, "%s%s" % ( ucsc_chrom_fa_filename, ext ) )
    
    fasta_base_filename = "%s.fa" % sequence_id
    fasta_filename = os.path.join( target_directory, fasta_base_filename )
    fasta_writer = open( fasta_filename, 'wb+' )
    
    tmp_extract_dir = os.path.join ( tmp_dir, 'extracted_fasta' )
    os.mkdir( tmp_extract_dir )
    
    tmp_fasta = open( ucsc_fasta_filename, 'wb+' )
    
    ftp.retrbinary( 'RETR %s' % ucsc_file_name, tmp_fasta.write )
    
    tmp_fasta.flush()
    tmp_fasta.seek( 0 )
    
    fasta_readers = get_stream_reader( tmp_fasta, tmp_extract_dir )
    
    data_table_entry = _stream_fasta_to_file( fasta_readers, target_directory, dbkey, sequence_id, sequence_name, params )
    _add_data_table_entry( data_manager_dict, data_table_entry )
    
    for fasta_reader in fasta_readers:
        fasta_reader.close()
    tmp_fasta.close()
    cleanup_before_exit( tmp_dir )

def download_from_ncbi( data_manager_dict, params, target_directory, dbkey, sequence_id, sequence_name ):
    NCBI_DOWNLOAD_URL = 'http://togows.dbcls.jp/entry/ncbi-nucleotide/%s.fasta' #FIXME: taken from dave's genome manager...why some japan site?
    
    requested_identifier = params['param_dict']['reference_source']['requested_identifier']
    url = NCBI_DOWNLOAD_URL % requested_identifier
    fasta_reader = urllib2.urlopen( url )
    
    data_table_entry = _stream_fasta_to_file( fasta_reader, target_directory, dbkey, sequence_id, sequence_name, params )
    _add_data_table_entry( data_manager_dict, data_table_entry )

def download_from_url( data_manager_dict, params, target_directory, dbkey, sequence_id, sequence_name ):
    #TODO: we should automatically do decompression here
    urls = filter( bool, map( lambda x: x.strip(), params['param_dict']['reference_source']['user_url'].split( '\n' ) ) )
    fasta_reader = [ urllib2.urlopen( url ) for url in urls ]
    
    data_table_entry = _stream_fasta_to_file( fasta_reader, target_directory, dbkey, sequence_id, sequence_name, params )
    _add_data_table_entry( data_manager_dict, data_table_entry )

def download_from_history( data_manager_dict, params, target_directory, dbkey, sequence_id, sequence_name ):
    #TODO: allow multiple FASTA input files
    input_filename = params['param_dict']['reference_source']['input_fasta']
    if isinstance( input_filename, list ):
        fasta_reader = [ open( filename, 'rb' ) for filename in input_filename ]
    else:
        fasta_reader = open( input_filename )
    
    data_table_entry = _stream_fasta_to_file( fasta_reader, target_directory, dbkey, sequence_id, sequence_name, params )
    _add_data_table_entry( data_manager_dict, data_table_entry )

def copy_from_directory( data_manager_dict, params, target_directory, dbkey, sequence_id, sequence_name ):
    input_filename = params['param_dict']['reference_source']['fasta_filename']
    create_symlink = params['param_dict']['reference_source']['create_symlink'] == 'create_symlink'
    if create_symlink:
        data_table_entry = _create_symlink( input_filename, target_directory, dbkey, sequence_id, sequence_name )
    else:
        if isinstance( input_filename, list ):
            fasta_reader = [ open( filename, 'rb' ) for filename in input_filename ]
        else:
            fasta_reader = open( input_filename )    
        data_table_entry = _stream_fasta_to_file( fasta_reader, target_directory, dbkey, sequence_id, sequence_name, params )
    _add_data_table_entry( data_manager_dict, data_table_entry )

def _add_data_table_entry( data_manager_dict, data_table_entry ):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables']['all_fasta'] = data_manager_dict['data_tables'].get( 'all_fasta', [] )
    data_manager_dict['data_tables']['all_fasta'].append( data_table_entry )
    return data_manager_dict

def _stream_fasta_to_file( fasta_stream, target_directory, dbkey, sequence_id, sequence_name, params, close_stream=True ):
    fasta_base_filename = "%s.fa" % sequence_id
    fasta_filename = os.path.join( target_directory, fasta_base_filename )
    fasta_writer = open( fasta_filename, 'wb+' )
    
    if isinstance( fasta_stream, list ) and len( fasta_stream ) == 1:
        fasta_stream = fasta_stream[0]
    
    if isinstance( fasta_stream, list ):
        last_char = None
        for fh in fasta_stream:
            if last_char not in [ None, '\n', '\r' ]:
                fasta_writer.write( '\n' )
            while True:
                data = fh.read( CHUNK_SIZE )
                if data:
                    fasta_writer.write( data )
                    last_char = data[-1]
                else:
                    break
            if close_stream:
                fh.close()
    else:
        while True:
            data = fasta_stream.read( CHUNK_SIZE )
            if data:
                fasta_writer.write( data )
            else:
                break
        if close_stream:
            fasta_stream.close()
    
    fasta_writer.close()
    
    sort_fasta( fasta_filename, params['param_dict']['sorting']['sort_selector'], params )
    
    return dict( value=sequence_id, dbkey=dbkey, name=sequence_name, path=fasta_base_filename )

def _create_symlink( input_filename, target_directory, dbkey, sequence_id, sequence_name ):
    fasta_base_filename = "%s.fa" % sequence_id
    fasta_filename = os.path.join( target_directory, fasta_base_filename )
    os.symlink( input_filename, fasta_filename )
    return dict( value=sequence_id, dbkey=dbkey, name=sequence_name, path=fasta_base_filename )




REFERENCE_SOURCE_TO_DOWNLOAD = dict( ucsc=download_from_ucsc, ncbi=download_from_ncbi, url=download_from_url, history=download_from_history, directory=copy_from_directory )

SORTING_METHODS = dict( as_is=_sort_fasta_as_is, lexicographical=_sort_fasta_lexicographical, gatk=_sort_fasta_gatk, custom=_sort_fasta_custom )

def main():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-d', '--dbkey_description', dest='dbkey_description', action='store', type="string", default=None, help='dbkey_description' )
    (options, args) = parser.parse_args()
    
    filename = args[0]
    
    params = loads( open( filename ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path']
    os.mkdir( target_directory )
    data_manager_dict = {}
    
    dbkey, sequence_id, sequence_name = get_dbkey_id_name( params, dbkey_description=options.dbkey_description ) 
    
    if dbkey in [ None, '', '?' ]:
        raise Exception( '"%s" is not a valid dbkey. You must specify a valid dbkey.' % ( dbkey ) )
    
    #Fetch the FASTA
    REFERENCE_SOURCE_TO_DOWNLOAD[ params['param_dict']['reference_source']['reference_source_selector'] ]( data_manager_dict, params, target_directory, dbkey, sequence_id, sequence_name )
    
    #save info to json file
    open( filename, 'wb' ).write( dumps( data_manager_dict ) )
        
if __name__ == "__main__": main()
