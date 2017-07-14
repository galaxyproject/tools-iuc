#!/usr/bin/env python

from __future__ import print_function
import optparse
import os.path
import sys
from query_db import get_connection, run_query, describe_tables


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
    parser.add_option('-o', '--output', dest='output', default=None,
                      help='Output file for query results')
    (options, args) = parser.parse_args()

    # determine output destination
    if options.output is not None:
        try:
            outputPath = os.path.abspath(options.output)
            outputFile = open(outputPath, 'w')
        except Exception as e:
            print("failed: %s" % e, file=sys.stderr)
            exit(3)
    else:
        outputFile = sys.stdout

    query = None
    if (options.query_file is not None):
        with open(options.query_file, 'r') as fh:
            query = ''
            for line in fh:
                query += line
    elif (options.query is not None):
        query = options.query

    if (query is None):
        try:
            describe_tables(get_connection(options.sqlitedb), outputFile)
        except Exception as exc:
            print("Error: %s" % exc, file=sys.stderr)
        exit(0)
    else:
        try:
            run_query(get_connection(options.sqlitedb), query, outputFile,
                      no_header=options.no_header)
        except Exception as exc:
            print("Error: %s" % exc, file=sys.stderr)
            exit(1)


if __name__ == "__main__":
    __main__()
