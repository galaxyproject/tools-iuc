#!/usr/bin/env python

import json
import optparse


def _add_data_table_entry( data_manager_dict, data_table_name, data_table_entry ):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables'][ data_table_name ] = data_manager_dict['data_tables'].get( data_table_name, [] )
    data_manager_dict['data_tables'][ data_table_name ].append( data_table_entry )
    return data_manager_dict


def build_rna_star_index(data_manager_dict, fasta_filename, target_directory,
                        dbkey, sequence_id, sequence_name, data_table_name,
                        sjdbOverhang, sjdbFileChrStartEnd,
                        sjdbGTFtagExonParentTranscript, sjdbGTFfeatureExon,
                        sjdbGTFchrPrefix):

    data_table_entry = dict( value=sequence_id, dbkey=dbkey, name=sequence_name, path=target_directory )
    data_manager_dict = _add_data_table_entry( data_manager_dict, data_table_name, data_table_entry )
    return data_manager_dict


def main():
    # Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '--fasta_filename', dest='fasta_filename', action='store', type="string", default=None, help='fasta_filename' )
    parser.add_option( '--fasta_dbkey', dest='fasta_dbkey', action='store', type="string", default=None, help='fasta_dbkey' )
    parser.add_option( '--fasta_description', dest='fasta_description', action='store', type="string", default=None, help='fasta_description' )
    parser.add_option( '--data_table_name', dest='data_table_name', action='store', type="string", default=None, help='data_table_name' )
    parser.add_option( '--sjdbGTFchrPrefix', type="string", default=None )
    parser.add_option( '--sjdbGTFfeatureExon', type="string", default=None )
    parser.add_option( '--sjdbGTFtagExonParentTranscript', type="string", default=None )
    parser.add_option( '--sjdbFileChrStartEnd', type="string", default=None )
    parser.add_option( '--sjdbOverhang', type="string", default='100' )
    (options, args) = parser.parse_args()
    filename = args[0]
    params = json.loads( open( filename ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path'].encode('ascii', 'replace')
    dbkey = options.fasta_dbkey
    if dbkey in [ None, '', '?' ]:
        raise Exception( '"%s" is not a valid dbkey. You must specify a valid dbkey.' % ( dbkey ) )

    sequence_id, sequence_name = options.fasta_dbkey, options.fasta_description

    data_manager_dict = build_rna_star_index(
        data_manager_dict={}, fasta_filename=options.fasta_filename,
        target_directory=target_directory, dbkey=dbkey, sequence_id=sequence_id,
        sequence_name=sequence_name, data_table_name=options.data_table_name,
        sjdbOverhang=options.sjdbOverhang,
        sjdbFileChrStartEnd=options.sjdbFileChrStartEnd,
        sjdbGTFtagExonParentTranscript=options.sjdbGTFtagExonParentTranscript,
        sjdbGTFfeatureExon=options.sjdbGTFfeatureExon,
        sjdbGTFchrPrefix=options.sjdbGTFchrPrefix)
    open( filename, 'wb' ).write( json.dumps( data_manager_dict ) )


if __name__ == "__main__":
    main()
