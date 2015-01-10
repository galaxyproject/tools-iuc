#!/usr/bin/env python
#Dan Blankenberg
# adapted from Dan's BWA one for rna star
# ross lazarus sept 2014
# fixed some stupid bugs January 2015


import sys
import os
import tempfile
import optparse
import subprocess

from galaxy.util.json import from_json_string, to_json_string


CHUNK_SIZE = 2**20
ONE_GB = 2**30

# from json import loads, dumps

DEFAULT_DATA_TABLE_NAME = "rnastar_index"


def get_id_name( params, dbkey, fasta_description=None):
    #TODO: ensure sequence_id is unique and does not already appear in location file
    sequence_id = params['param_dict']['sequence_id']
    if not sequence_id:
        sequence_id = dbkey
    
    sequence_name = params['param_dict']['sequence_name']
    if not sequence_name:
        sequence_name = fasta_description
        if not sequence_name:
            sequence_name = dbkey
    return sequence_id, sequence_name


def _run_command( command, target_directory ):
    tmp_stderr = tempfile.NamedTemporaryFile( prefix = "tmp-data-manager-rnastar_index_builder-stderr" )
    return_code = subprocess.call( command, shell=False, stdout=sys.stdout, stderr=tmp_stderr.fileno() )
    if return_code:
        tmp_stderr.flush()
        tmp_stderr.seek( 0 )
        sys.stderr.write( "### star reports an error:\n" )
        while True:
            chunk = tmp_stderr.read( CHUNK_SIZE )
            if not chunk:
                break
            sys.stderr.write( chunk )
        sys.exit( return_code )
    tmp_stderr.close()
    return return_code


def build_rnastar_index( data_manager_dict, fasta_filename,  target_directory, 
    dbkey, sequence_id, sequence_name, data_table_name,
    sjdbOverhang,sjdbGTFfile, sjdbFileChrStartEnd,sjdbGTFtagExonParentTranscript,sjdbGTFfeatureExon,sjdbGTFchrPrefix,n_threads):
    #TODO: allow multiple FASTA input files
    fasta_base_name = os.path.basename( fasta_filename )
    sym_linked_fasta_filename = os.path.join( target_directory, fasta_base_name )
    os.symlink( fasta_filename, sym_linked_fasta_filename )
    print >> sys.stdout,'made',sym_linked_fasta_filename
    cl = ['STAR','--runMode', 'genomeGenerate', '--genomeFastaFiles', sym_linked_fasta_filename, '--genomeDir', target_directory, '--runThreadN',n_threads ]
    
    if sjdbGTFfile:
         cl += [ '--sjdbGTFfeatureExon', sjdbGTFfeatureExon,'--sjdbGTFtagExonParentTranscript',sjdbGTFtagExonParentTranscript]
         if (sjdbGTFchrPrefix > ''):
             cl += ['--sjdbGTFchrPrefix', sjdbGTFchrPrefix]
         cl += ['--sjdbOverhang', sjdbOverhang, '--sjdbGTFfile', sjdbGTFfile]
    elif sjdbFileChrStartEnd:
        cl += ['--sjdbFileChrStartEnd', sjdbFileChrStartEnd, '--sjdbOverhang', sjdbOverhang]

    tmp_stderr = tempfile.NamedTemporaryFile( prefix = "tmp-data-manager-rnastar-index-builder-stderr" )
    proc = subprocess.Popen( args=cl, shell=False, cwd=target_directory, stderr=tmp_stderr.fileno() )
    return_code = proc.wait()
    if return_code:
        tmp_stderr.flush()
        tmp_stderr.seek(0)
        print >> sys.stderr, "Error building index:"
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
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '--fasta_filename', dest='fasta_filename', action='store', type="string", default=None, help='fasta_filename' )
    parser.add_option( '--fasta_dbkey', dest='fasta_dbkey', action='store', type="string", default=None, help='fasta_dbkey' )
    parser.add_option( '--fasta_description', dest='fasta_description', action='store', type="string", default=None, help='fasta_description' )
    parser.add_option( '--data_table_name', dest='data_table_name', action='store', type="string", default=None, help='data_table_name' )
    parser.add_option( '--sjdbGTFfile', type="string", default=None )
    parser.add_option( '--sjdbGTFchrPrefix', type="string", default=None )
    parser.add_option( '--sjdbGTFfeatureExon', type="string", default=None )
    parser.add_option( '--sjdbGTFtagExonParentTranscript', type="string", default=None )
    parser.add_option( '--sjdbFileChrStartEnd', type="string", default=None )
    parser.add_option( '--sjdbOverhang', type="string", default='100' )
    parser.add_option( '--runThreadN', type="string", default='4' )
    (options, args) = parser.parse_args()
    filename = args[0]
    params = from_json_string( open( filename ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path'].encode('ascii','replace')
    data_manager_dict = {}
    dbkey = options.fasta_dbkey
    if dbkey in [ None, '', '?' ]:
        raise Exception( '"%s" is not a valid dbkey. You must specify a valid dbkey.' % ( dbkey ) )
    
    sequence_id, sequence_name = get_id_name( params, dbkey=dbkey, fasta_description=options.fasta_description )
    
    try:
         os.mkdir( target_directory )
    except OSError:
        pass
    #build the index
    build_rnastar_index( data_manager_dict = data_manager_dict, fasta_filename = options.fasta_filename, target_directory = target_directory,
      dbkey = dbkey, sequence_id = sequence_id, sequence_name = sequence_name, data_table_name=options.data_table_name,
      sjdbOverhang=options.sjdbOverhang,sjdbGTFfile=options.sjdbGTFfile,
      sjdbFileChrStartEnd=options.sjdbFileChrStartEnd,sjdbGTFtagExonParentTranscript=options.sjdbGTFtagExonParentTranscript,
      sjdbGTFfeatureExon=options.sjdbGTFfeatureExon,sjdbGTFchrPrefix=options.sjdbGTFchrPrefix,
      n_threads=options.runThreadN )
    open( filename, 'wb' ).write( to_json_string( data_manager_dict ) )
        
if __name__ == "__main__": main()

