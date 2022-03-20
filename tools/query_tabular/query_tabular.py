#!/usr/bin/env python

from __future__ import print_function

import json
import optparse
import os.path
import sys

from load_db import create_table
from query_db import describe_tables, get_connection, run_query


"""
JSON config:
{ tables : [
    { file_path : '/home/galaxy/dataset_101.dat',
            table_name : 't1',
            column_names : ['c1','c2','c3'],
            pkey_autoincr : 'id'
            comment_lines : 1
            unique: ['c1'],
            index: ['c2', 'c3']
    },
    { file_path : '/home/galaxy/dataset_102.dat',
            table_name : 'gff',
            column_names : ['seqname',,'date','start','end']
            comment_lines : 1
            load_named_columns : True
            filters : [{'filter': 'regex', 'pattern': '#peptide',
                        'action': 'exclude_match'},
                       {'filter': 'replace', 'column': 3,
                        'replace': 'gi[|]', 'pattern': ''}]
    },
    { file_path : '/home/galaxy/dataset_103.dat',
            table_name : 'test',
            column_names : ['c1', 'c2', 'c3']
    }
    ]
}
"""


def __main__():
    # Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option('-s', '--sqlitedb', dest='sqlitedb', default=None,
                      help='The SQLite Database')
    parser.add_option('-j', '--jsonfile', dest='jsonfile', default=None,
                      help='JSON dict of table specifications')
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
    parser.add_option('-d', '--debug', dest='debug', default=False,
                      action='store_true',
                      help='Output info to stderr')
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

    def _create_table(ti, table):
        path = table['file_path']
        table_name =\
            table['table_name'] if 'table_name' in table else 't%d' % (ti + 1)
        comment_lines =\
            table['comment_lines'] if 'comment_lines' in table else 0
        comment_char =\
            table['comment_char'] if 'comment_char' in table else None
        column_names =\
            table['column_names'] if 'column_names' in table else None
        firstlinenames =\
            table['firstlinenames'] if 'firstlinenames' in table else False
        if column_names:
            load_named_columns =\
                table['load_named_columns']\
                if 'load_named_columns' in table else False
        else:
            load_named_columns = False
        unique_indexes = table['unique'] if 'unique' in table else []
        indexes = table['index'] if 'index' in table else []
        filters = table['filters'] if 'filters' in table else None
        pkey_autoincr = \
            table['pkey_autoincr'] if 'pkey_autoincr' in table else None
        create_table(get_connection(options.sqlitedb), path, table_name,
                     pkey_autoincr=pkey_autoincr,
                     firstlinenames=firstlinenames,
                     column_names=column_names,
                     skip=comment_lines,
                     comment_char=comment_char,
                     load_named_columns=load_named_columns,
                     filters=filters,
                     unique_indexes=unique_indexes,
                     indexes=indexes)

    if options.jsonfile:
        try:
            with open(options.jsonfile) as fh:
                tdef = json.load(fh)
                if options.debug:
                    print('JSON: %s' % tdef, file=sys.stderr)
                if 'tables' in tdef:
                    for ti, table in enumerate(tdef['tables']):
                        _create_table(ti, table)
                if 'sql_stmts' in tdef:
                    for si, stmt in enumerate(tdef['sql_stmts']):
                        rowcount = run_query(get_connection(options.sqlitedb), stmt, None)
                        if options.debug:
                            print('\nDB modification: %s  \nrowcount: %s' %
                                  (stmt, rowcount), file=sys.stderr)
                if 'queries' in tdef:
                    for qi, qstmt in enumerate(tdef['queries']):
                        if 'header' in qstmt:
                            no_header = False
                            comment_char = qstmt['header']
                        else:
                            no_header = True
                            comment_char = None
                        with open(qstmt['result_file'], 'w') as fh:
                            query = qstmt['query']
                            rowcount = run_query(get_connection(options.sqlitedb),
                                                 query,
                                                 fh,
                                                 no_header=no_header,
                                                 comment_char=comment_char)
                        if options.debug:
                            print('\nSQL: %s  \nrowcount: %s' %
                                  (query, rowcount), file=sys.stderr)
        except Exception as e:
            exit('Error: %s' % (e))

    query = None
    if options.query_file is not None:
        with open(options.query_file, 'r') as fh:
            query = ''
            for line in fh:
                query += line
    elif options.query is not None:
        query = options.query

    if query is None:
        try:
            describe_tables(get_connection(options.sqlitedb), outputFile)
        except Exception as e:
            exit('Error: %s' % (e))
    else:
        try:
            rowcount = run_query(get_connection(options.sqlitedb),
                                 query, outputFile,
                                 no_header=options.no_header,
                                 comment_char=options.comment_char)
            if options.debug:
                print('\nSQL: %s  \nrowcount: %s' %
                      (query, rowcount), file=sys.stderr)
        except Exception as e:
            exit('Error: %s' % (e))


if __name__ == "__main__":
    __main__()
