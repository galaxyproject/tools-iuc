#!/usr/bin/env python

from __future__ import print_function

import optparse
import os.path
import sys

from query_db import describe_tables, get_connection, run_query


def __main__():
    # Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option('-s', '--sqlitedb', dest='sqlitedb', default=None,
                      help='The SQLite Database')
    parser.add_option('-q', '--query', dest='query', default=None,
                      help='SQL query')
    parser.add_option('-Q', '--query_file', dest='query_file', default=None,
                      help='SQL query file')
    parser.add_option('-n', '--no_header', dest='no_header', default=False,
                      action='store_true',
                      help='Include a column headers line')
    parser.add_option('-c', '--comment_char', dest='comment_char', default='',
                      help='comment character to prefix column header line')
    parser.add_option('-o', '--output', dest='output', default=None,
                      help='Output file for query results')
    (options, args) = parser.parse_args()

    # determine output destination
    if options.output is not None:
        try:
            outputPath = os.path.abspath(options.output)
            outputFile = open(outputPath, 'w')
        except Exception as e:
            exit('Error: %s' % (e))
    else:
        outputFile = sys.stdout

    query = None
    if options.query_file is not None:
        with open(options.query_file, 'r') as fh:
            query = fh.read()
    elif options.query is not None:
        query = options.query

    if query is None:
        try:
            describe_tables(get_connection(options.sqlitedb), outputFile)
        except Exception as e:
            exit('Error: %s' % (e))
        exit(0)
    else:
        try:
            run_query(get_connection(options.sqlitedb), query, outputFile,
                      no_header=options.no_header,
                      comment_char=options.comment_char)
        except Exception as e:
            exit('Error: %s' % (e))


if __name__ == "__main__":
    __main__()
