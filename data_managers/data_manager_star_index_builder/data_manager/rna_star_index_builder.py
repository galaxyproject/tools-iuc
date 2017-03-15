#!/usr/bin/env python

import json
import optparse




def main():
    # Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '--fasta_filename', dest='fasta_filename', action='store', type="string", default=None, help='fasta_filename' )
    parser.add_option( '--fasta_dbkey', dest='fasta_dbkey', action='store', type="string", default=None, help='fasta_dbkey' )
    parser.add_option( '--fasta_description', dest='fasta_description', action='store', type="string", default=None, help='fasta_description' )
    parser.add_option( '--data_table_name', dest='data_table_name', action='store', type="string", default=None, help='data_table_name' )
    parser.add_option( '--subdir', dest='subdir', action='store', type="string", default=None, help='subdir' )
    (options, args) = parser.parse_args()
    filename = args[0]
    params = json.loads( open( filename ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path'].encode('ascii', 'replace')
    dbkey = options.fasta_dbkey
    if dbkey in [ None, '', '?' ]:
        raise Exception( '"%s" is not a valid dbkey. You must specify a valid dbkey.' % ( dbkey ) )



    data_manager_dict = {'data_tables': {options.data_table_name: [dict( value=dbkey, dbkey=dbkey, name=options.fasta_description, path=options.subdir )]}}

    open( filename, 'wb' ).write( json.dumps( data_manager_dict ) )


if __name__ == "__main__":
    main()
