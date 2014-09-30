#!/usr/bin/env python
#Dan Blankenberg
# adapted from Dan's BWA one for rna star
# ross lazarus sept 2014
#
import sys
import os
import tempfile
import optparse
import subprocess

from galaxy.util.json import from_json_string, to_json_string

DEFAULT_DATA_TABLE_NAME = "rnastar_indexes"

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

def build_rnastar_index( data_manager_dict, fasta_filename, target_directory, dbkey, sequence_id, sequence_name, data_table_name,
    sjdbOverhang,sjdbGTFfile, sjdbFileChrStartEnd,sjdbGTFtagExonParentTranscript,sjdbGTFfeatureExon,sjdbGTFchrPrefix,n_threads):
    #TODO: allow multiple FASTA input files
    #tmp_dir = tempfile.mkdtemp( prefix='tmp-data-manager-bwa-index-builder-' )
    fasta_base_name = os.path.split( fasta_filename )[-1]
    sym_linked_fasta_filename = os.path.join( target_directory, fasta_base_name )
    os.symlink( fasta_filename, sym_linked_fasta_filename )
    pdict={'target_directory':target_directory,'n_threads':n_threads, 'sjdbFileChrStartEnd':sjdbFileChrStartEnd,
           'sjdbGTFtagExonParentTranscript':sjdbGTFtagExonParentTranscript, 'sjdbGTFfeatureExon':sjdbGTFfeatureExon,
           'sjdbGTFchrPrefix':sjdbGTFchrPrefix,'sjdbOverhang':sjdbOverhang, 'sjdbGTFfile':sjdbGTFfile,
           'sym_linked_fasta_filename':sym_linked_fasta_filename}
    
    cl = 'STAR --runMode genomeGenerate --genomeFastaFiles %(sym_linked_fasta_filename)s --genomeDir %(target_directory)s --runThreadN %(n_threads)s' % pdict
    if sjdbGTFfile:
         cl += '''--sjdbGTFchrPrefix %(sjdbGTFchrPrefix)s --sjdbGTFfeatureExon %(sjdbGTFfeatureExon)s --sjdbOverhang %(sjdbOverhang)s
   --sjdbGTFfile %(sjdbGTFfile)s --sjdbGTFtagExonParentTranscript %(sjdbGTFtagExonParentTranscript)s''' %  pdict
    elif sjdbFileChrStartEnd:
        cl += '--sjdbFileChrStartEnd %(sjdbFileChrStartEnd)s --sjdbOverhang %(sjdbOverhangs)s' % pdict
    tmp_stderr = tempfile.NamedTemporaryFile( prefix = "tmp-data-manager-rnastar-index-builder-stderr" )
    args = cl.split(' ')
    proc = subprocess.Popen( args=args, shell=False, cwd=target_directory, stderr=tmp_stderr.fileno() )
    return_code = proc.wait()
    if return_code:
        tmp_stderr.flush()
        tmp_stderr.seek(0)
        print >> sys.stderr, "Error building index: retcode=",retcode
        while True:
            chunk = tmp_stderr.read( CHUNK_SIZE )
            if not chunk:
                break
            sys.stderr.write( chunk )
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
    parser.add_option( '-f', '--fasta_filename', dest='fasta_filename', action='store', type="string", default=None, help='fasta_filename' )
    parser.add_option( '-d', '--fasta_dbkey', dest='fasta_dbkey', action='store', type="string", default=None, help='fasta_dbkey' )
    parser.add_option( '-t', '--fasta_description', dest='fasta_description', action='store', type="string", default=None, help='fasta_description' )
    parser.add_option( '-n', '--data_table_name', dest='data_table_name', action='store', type="string", default=None, help='data_table_name' )
    parser.add_option( '--out_file', default=None)
    parser.add_option( '--out_index_path', default=None)
    parser.add_option( '--sjdbGTFfile', type="string", default=None )
    parser.add_option( '--sjdbGTFchrPrefix', type="string", default=None )
    parser.add_option( '--sjdbGTFfeatureExon', type="string", default=None )
    parser.add_option( '--sjdbGTFtagExonParentTranscript', type="string", default=None )
    parser.add_option( '--sjdbFileChrStartEnd', type="string", default=None )
    parser.add_option( '--sjdbOverhang', type="int", default=100 )
    parser.add_option( '--runThreadN', type="int", default=4 )
    (options, args) = parser.parse_args()
    
    filename = options.out_file
    params = from_json_string( open( filename ).read() )
    target_directory = options.out_index_path
    os.mkdir( target_directory )
    data_manager_dict = {}
    
    dbkey = options.fasta_dbkey
    
    if dbkey in [ None, '', '?' ]:
        raise Exception( '"%s" is not a valid dbkey. You must specify a valid dbkey.' % ( dbkey ) )
    
    sequence_id, sequence_name = get_id_name( params, dbkey=dbkey, fasta_description=options.fasta_description )
    
    #build the index
    build_rnastar_index( data_manager_dict, options.fasta_filename, target_directory, dbkey, sequence_id, sequence_name, data_table_name=options.data_table_name,
      sjdbOverhang=options.sjdbOverhang,sjdbGTFfile=options.sjdbGTFfile,
      sjdbFileChrStartEnd=options.sjdbFileChrStartEnd,sjdbGTFtagExonParentTranscript=options.sjdbGTFtagExonParentTranscript,
      sjdbGTFfeatureExon=options.sjdbGTFfeatureExon,sjdbGTFchrPrefix=options.sjdbGTFchrPrefix,
      n_threads=options.runThreadN )

    
    #save info to json file
    open( filename, 'wb' ).write( to_json_string( data_manager_dict ) )
        
if __name__ == "__main__": main()
