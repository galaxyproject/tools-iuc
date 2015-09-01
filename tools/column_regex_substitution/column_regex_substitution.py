#!/usr/bin/env python
#Dan Blankenberg

"""
A script for using regex substitutions on columns.
"""

import optparse
import re
import sys
import string

VERSION = "0.0.1"

COLUMN_STRIP_VALUES = "".join( set( string.printable ) - set( string.digits ) - set(',') )

def get_provided_columns( provided_value, column_offset ):
    try:
        rval = sorted( map( lambda x: int( x.strip( COLUMN_STRIP_VALUES ) ) + column_offset, provided_value.split( ',' ) ) )
    except:
        rval = None
    if rval:
        return rval
    return None


def __main__():
    parser = optparse.OptionParser()
    parser.add_option('--pattern', action='store', default=None,
                      help='pattern string')
    parser.add_option('--replacement', action='store', default=None,
                      help='replacement string')
    parser.add_option('--input', action='store', default=None,
                      help='Filename of input file')
    parser.add_option('--output', action='store', default=None,
                      help='Filename of output file')
    parser.add_option('--delimiter', action='store', default=None,
                      help='column delimiter')
    parser.add_option('--columns', action='store', default=None,
                      help='columns to operate on')
    parser.add_option('--column_offset', action='store', default=0,
                      help='offset to apply to columns index to force to zero-based')
    parser.add_option('--skip', action='store', default=0,
                      help='Number of lines to skip')
    parser.add_option('--version', action='store_true', default=False,
                      help='Show version')

    (options, args) = parser.parse_args()

    if options.version:
        print "blankenberg_python_regex_substitution %s" % ( VERSION )
        sys.exit(0)

    if None in [ options.pattern, options.replacement, options.output ]:
        parser.print_help()
        sys.exit(1)

    pattern = options.pattern
    replacement = options.replacement
    column_offset = int( options.column_offset )
    print "Pattern: %s\nReplacement: %s" % ( repr( pattern ), repr( replacement ) )
    pattern = re.compile( pattern )
    provided_columns = get_provided_columns( options.columns, column_offset )
    if provided_columns:
        column_str = ", ".join( map( lambda x: str( x - column_offset ), provided_columns ) )
    else:
        column_str = 'all'
    print "With delimiter %s, on columns: %s" % ( repr( options.delimiter ), column_str )
    if options.delimiter is None:
        split_func = lambda x: [ x.rstrip( '\n\r' ) ]
        join_char = ""
    else:
        split_func = lambda x: x.rstrip( '\n\r' ).split( options.delimiter )
        join_char = options.delimiter
    with open( options.input, 'rb' ) as fin:
        with open( options.output, 'w') as fout:
            for i, line in enumerate( fin ):
                if i < options.skip:
                    continue
                line = split_func( line )
                field_count = len( line )
                if provided_columns:
                    columns = provided_columns
                else:
                    columns = range( field_count )
                for j in columns:
                    if j >= field_count:
                        break
                    line[ j ] = re.sub( pattern, replacement, line[ j ] )
                fout.write( "%s\n" % ( join_char.join( line ) ) )

if __name__ == "__main__":
    __main__()
