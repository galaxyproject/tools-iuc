#!/usr/bin/env python
"""
This tool takes a tab-delimited textfile as input and creates another column in
the file which is the result of a computation performed on every row in the
original file. The tool will skip over invalid lines within the file,
informing the user about the number of lines skipped.
"""

import argparse
import json
import re
import sys
# functions that may be used in the compute expression
from math import (  # noqa: F401
    ceil,
    exp,
    floor,
    log,
    log10,
    sqrt
)

from numpy import format_float_positional  # noqa: F401

parser = argparse.ArgumentParser()
parser.add_argument('input', type=argparse.FileType('r'), help="input file")
parser.add_argument('output', type=argparse.FileType('wt'), help="output file")
parser.add_argument('cond', nargs='?', type=str, help="expression")
parser.add_argument('columns', nargs='?', type=int, help="number of columns")
parser.add_argument('column_types', nargs='?', type=str, help="comma separated list of column types")
parser.add_argument('--round', action="store_true",
                    help="round result")
parser.add_argument('--avoid_scientific_notation', action="store_true",
                    help="avoid scientific notation")
parser.add_argument('--header_new_column_name', default=None, type=str,
                    help="First line of input is a header line with column "
                         "names and this should become the name of the new "
                         "column")
parser.add_argument('--load_json', default=None, type=argparse.FileType('r'),
                    help="overwrite parsed arguments from json file")
args = parser.parse_args()

argparse_dict = vars(args)
if args.load_json:
    json_dict = json.load(args.load_json)
    argparse_dict.update(json_dict)

fh = argparse_dict['input']
out = argparse_dict['output']
expr = argparse_dict['cond']
round_result = argparse_dict['round']
avoid_scientific_notation = argparse_dict['avoid_scientific_notation']

if argparse_dict['header_new_column_name'] is not None:
    header_line = fh.readline().strip('\n')
    out.write(
        '{0}\t{1}\n'.format(
            header_line, argparse_dict['header_new_column_name']
        )
    )
try:
    in_columns = int(argparse_dict['columns'])
    if in_columns < 1:
        # To be considered tabular, data must have at least one column.
        raise ValueError
except Exception:
    if not fh.readline():
        # empty file content is ok and should produce empty output
        out.close()
        sys.exit()
    sys.exit("Missing or invalid 'columns' metadata value, click the pencil icon in the history item and select the Auto-detect option to correct it.  This tool can only be used with tab-delimited data.")
try:
    in_column_types = argparse_dict['column_types'].split(',')
except Exception:
    sys.exit("Missing or invalid 'column_types' metadata value, click the pencil icon in the history item and select the Auto-detect option to correct it.  This tool can only be used with tab-delimited data.")
if len(in_column_types) != in_columns:
    sys.exit("The 'columns' metadata setting does not conform to the 'column_types' metadata setting, click the pencil icon in the history item and select the Auto-detect option to correct it.  This tool can only be used with tab-delimited data.")

operators = 'is|not|or|and'
builtin_and_math_functions = 'abs|all|any|bin|chr|cmp|complex|divmod|float|bool|hex|int|len|long|max|min|oct|ord|pow|range|reversed|round|sorted|str|sum|type|unichr|unicode|log|log10|exp|sqrt|ceil|floor'
string_and_list_methods = [name for name in dir('') + dir([]) if not name.startswith('_')]
whitelist = r"^([c0-9\+\-\*\/\(\)\.\'\"><=,:! ]|%s|%s|%s)*$" % (operators, builtin_and_math_functions, '|'.join(string_and_list_methods))
if not re.compile(whitelist).match(expr):
    sys.exit("Invalid expression")
if avoid_scientific_notation:
    expr = "format_float_positional(%s)" % expr

# Prepare the column variable names and wrappers for column data types
cols, type_casts = [], []
for col in range(1, in_columns + 1):
    col_name = "c%d" % col
    cols.append(col_name)
    col_type = in_column_types[col - 1].strip()
    if not round_result and col_type == 'int':
        col_type = 'float'
    type_cast = "%s(%s)" % (col_type, col_name)
    type_casts.append(type_cast)

col_str = ', '.join(cols)    # 'c1, c2, c3, c4'
type_cast_str = ', '.join(type_casts)  # 'str(c1), int(c2), int(c3), str(c4)'
assign = "%s = line.split('\\t')" % col_str
if len(cols) == 1:
    # Single column, unpacking by assignment won't work
    assign += '[0]'
wrap = "%s = %s" % (col_str, type_cast_str)
skipped_lines = 0
first_invalid_line = 0
invalid_line = None
lines_kept = 0
total_lines = 0

# Read input file, skipping invalid lines, and perform computation that will result in a new column
code = '''
for i, line in enumerate(fh):
    total_lines += 1
    line = line.rstrip('\\r\\n')
    if not line or line.startswith('#'):
        skipped_lines += 1
        if not invalid_line:
            first_invalid_line = i + 1
            invalid_line = line
        continue
    try:
        %s
        %s
        new_val = %s
        if round_result:
            new_val = int(round(new_val))
        new_line = line + '\\t' + str(new_val) + "\\n"
        out.write(new_line)
        lines_kept += 1
    except Exception:
        skipped_lines += 1
        if not invalid_line:
            first_invalid_line = i + 1
            invalid_line = line
fh.close()
''' % (assign, wrap, expr)

valid_expr = True
try:
    exec(code)
except Exception as e:
    out.close()
    if str(e).startswith('invalid syntax'):
        valid_expr = False
        sys.exit('Expression "%s" likely invalid. See tool tips, syntax and examples.' % expr)
    else:
        sys.exit(str(e))

if valid_expr:
    out.close()
    valid_lines = total_lines - skipped_lines
    print('Creating column %d with expression %s' % (in_columns + 1, expr))
    if valid_lines > 0:
        print('kept %4.2f%% of %d lines.' % (100.0 * lines_kept / valid_lines,
                                             total_lines))
    else:
        print('Possible invalid expression "%s" or non-existent column referenced. See tool tips, syntax and examples.' % expr)
    if skipped_lines > 0:
        print('Skipped %d invalid lines starting at line #%d: "%s"' %
              (skipped_lines, first_invalid_line, invalid_line))
