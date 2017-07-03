#!/usr/bin/env python

import json
import optparse


def main():
    parser = optparse.OptionParser()
    parser.add_option( '--config-file', dest='config_file', action='store', type="string")
    parser.add_option( '--value', dest='value', action='store', type="string" )
    parser.add_option( '--dbkey', dest='dbkey', action='store', type="string" )
    parser.add_option( '--name', dest='name', action='store', type="string" )
    parser.add_option( '--subdir', dest='subdir', action='store', type="string" )
    parser.add_option( '--data-table', dest='data_table', action='store', type="string" )
    parser.add_option( '--withGTF', dest='withGTF', action='store_true' )
    (options, args) = parser.parse_args()

    if options.dbkey in [ None, '', '?' ]:
        raise Exception( '"%s" is not a valid dbkey. You must specify a valid dbkey.' % ( options.dbkey ) )

    withGTF = "0"
    if options.withGTF:
        withGTF = "1"

    data_manager_dict = {'data_tables': {options.data_table: [dict({"value": options.value, "dbkey": options.dbkey, "name": options.name, "path": options.subdir, "with-gtf": withGTF} )]}}
    open( options.config_file, 'wb' ).write( json.dumps( data_manager_dict ) )


if __name__ == "__main__":
    main()
