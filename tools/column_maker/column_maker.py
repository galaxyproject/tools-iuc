#!/usr/bin/env python
"""
This tool takes a tab-delimited textfile as input and creates new columns in
the file which are the result of a computation performed on every row in the
original file. The tool will skip over invalid lines within the file,
informing the user about the number of lines skipped. It does not change the
formatting of any original, retained columns.
"""

import argparse
import enum
import re
import sys
# Functions that may be used in the compute expression
from math import (  # noqa: F401
    ceil,
    exp,
    floor,
    log,
    log10,
    sqrt,
)

from numpy import format_float_positional


class Mode(enum.Enum):
    APPEND = ''
    INSERT = 'I'
    REPLACE = 'R'


parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, help='input file')
parser.add_argument('output', type=str, help='output file')
parser.add_argument(
    '-t', '--column-types', nargs='?', required=True,
    help='A comma-separated list of column types in the input file'
)
parser.add_argument(
    '--avoid-scientific-notation', action='store_true',
    help='avoid scientific notation'
)
parser.add_argument(
    '--header', action='store_true',
    help='The input has a header line with column names. '
         'Expressions must specify names of newly calculated columns.'
)
parser.add_argument(
    '--fail-on-non-existent-columns', action='store_true',
    help='If an expression references a column number that is not existent '
         'when the expression gets computed, the default behavior is to treat '
         'this as a case of rows for which the expression cannot be computed. '
         'The behavior of the tool will then depend on which of the '
         'non-computable switches is in effect. With this flag, in contrast, '
         'the tool will fail directly upon encountering a non-existing column.'
)
non_computable = parser.add_mutually_exclusive_group()
non_computable.add_argument('--skip-non-computable', action='store_true')
non_computable.add_argument('--keep-non-computable', action='store_true')
non_computable.add_argument('--fail-on-non-computable', action='store_true')
non_computable.add_argument('--non-computable-blank', action='store_true')
non_computable.add_argument('--non-computable-default')

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-e', '--expr', nargs='*', type=str, help='Expression(s)')
group.add_argument(
    '-f', '--file', type=str,
    help='File to read expressions from (mutually exclusive with -e)'
)
args = parser.parse_args()

if not args.column_types:
    with open(args.input) as fh:
        if not fh.readline():
            # Generally, the input must have at least one column to be
            # considered tabular, but empty files are ok and should produce
            # empty output.
            with open(args.output, 'w') as out:
                pass
            sys.exit()
        sys.exit(
            "Missing column types. "
            "In Galaxy, click the pencil icon on the history item and "
            "select the Auto-detect option to correct it.  "
            "This tool can only be used with tab-delimited data."
        )

in_column_types = [t.strip() for t in args.column_types.split(',')]
in_columns = len(in_column_types)

# Prepare initial column variable names and type cast representations
# for column data types
cols, type_casts = [], []
for n, col_type in enumerate(in_column_types, start=1):
    col_name = "c%d" % n
    cols.append(col_name)
    type_cast = "%s(%s)" % (col_type, col_name)
    type_casts.append(type_cast)

col_str = ', '.join(cols)    # 'c1, c2, c3, c4'
type_cast_str = ', '.join(type_casts)  # 'str(c1), int(c2), int(c3), str(c4)'
# Define lambda for type-casting of original row fields
try:
    cast_types = eval('lambda %s: [%s]' % (col_str, type_cast_str))
except Exception as e:
    sys.exit(
        'While parsing column types, the following problem occured: "%s"'
        % e
    )

# Get and parse expressions
if args.file:
    expr = []
    with open(args.file) as i:
        for line in i:
            line = line.strip()
            if line:
                expr.append(line)
else:
    expr = args.expr

# each expr must be a full data row manipulation instruction of the form:
# EXPR;[COL_ADD_SPEC];[COL_NAME]
# where EXPR is the actual expression to compute on the row,
# COL_ADD_SPEC consists of a column index and a mode identifier for how the
# new column should be added.
# Examples: 3I (insert new col before current column 3),
# 2R (replace current column 2 with new column);
# a missing COL_ADD_SPEC is interpreted as mode A (append new column at the
# end of the row).
# COL_NAME is required with the --header option and specifies the name of the
# new column; without --header, any COL_NAME gets ignored.
operators = 'is|not|or|and'
builtin_and_math_functions = 'abs|all|any|bin|chr|cmp|complex|divmod|float|bool|hex|int|len|long|max|min|oct|ord|pow|range|reversed|round|sorted|str|sum|type|unichr|unicode|log|log10|exp|sqrt|ceil|floor'
imported_numpy_function = 'format_float_positional'
string_and_list_methods = [name for name in dir('') + dir([]) if not name.startswith('_')]
whitelist = r"^([c0-9\+\-\*\/\(\)\.\'\"><=,:! ]|%s|%s|%s|%s)*$" % (
    operators,
    builtin_and_math_functions,
    imported_numpy_function,
    '|'.join(string_and_list_methods)
)
valid_pat = re.compile(whitelist)
ops = []
num_cols = in_columns
for e in expr:
    parts = e.split(';')
    if len(parts) != 3:
        sys.exit(
            'Invalid compute instruction: "%s".  '
            'Required format: EXPR;[COL_ADD_SPEC];[COL_NAME]' % e
        )
    if not valid_pat.match(parts[0]):
        sys.exit('Invalid expression: "%s"' % parts[0])
    try:
        exp_lambda = eval('lambda %s: %s' % (col_str, parts[0]))
    except Exception as e:
        if str(e).startswith('invalid syntax'):
            sys.exit(
                'Expression "%s" caused a syntax error during parsing.'
                % parts[0]
            )
        else:
            sys.exit(
                'While parsing expression "%s" the following problem occured: '
                '"%s"' % (parts[0], str(e))
            )
    try:
        new_col_idx = int(parts[1][:-1] or '0') - 1
    except ValueError:
        sys.exit(
            'COL_ADD_SPECS need to start with a (1-based) column index. '
            'Could not parse a column index from "%s"' % parts[1]
        )
    try:
        mode = Mode(parts[1][-1:])
    except ValueError:
        sys.exit(
            'COL_ADD_SPECS need to end in a single-character mode identifier '
            '("A", "I", or "R").  '
            'Could not parse a valid identifier from "%s"' % parts[1]
        )
    if mode is Mode.APPEND and new_col_idx != -1:
        sys.exit(
            '"%s" combines Append mode with a column index.  '
            'Please specify mode "A" without a leading column index.'
            % parts[1]
        )
    if mode is Mode.REPLACE:
        if new_col_idx < 0 or new_col_idx >= num_cols:
            sys.exit(
                'Cannot replace the contents of column %d as specified by '
                'instruction "%s".  No such column at this point of the '
                'computation' % (new_col_idx + 1, e)
            )
    new_col_name = parts[2] if len(parts) == 3 else None
    if not new_col_name and args.header:
        sys.exit(
            'A name is required for any new columns when using an existing '
            'header line (--header option), but found none in instruction: '
            '"%s"' % e
        )
    ops.append([exp_lambda, new_col_idx, mode, new_col_name, parts[0]])
    if mode is Mode.APPEND or mode is Mode.INSERT:
        # If the current expression results in an additional column,
        # we need to handle the new field in subsequent lambda functions.
        num_cols += 1
        col_str += ', c%d' % num_cols


# ready to start parsing the input file
print(
    'Computing %d columns with instructions %s' % (num_cols - in_columns, expr)
)
skipped_lines = 0
first_invalid_line = 0
invalid_line = None
lines_computed = 0
total_lines = 0
non_existent_col_pat = re.compile(r"name 'c\d+' is not defined")

with open(args.input) as fh, open(args.output, 'w') as out:
    if args.header:
        # compute new header line from original
        header_cols = fh.readline().strip('\n').split('\t')
        for _, col_idx, mode, col_name, _ in ops:
            if mode is Mode.INSERT:
                header_cols.insert(col_idx, col_name)
            elif mode is Mode.REPLACE:
                header_cols[col_idx] = col_name
            else:
                header_cols.append(col_name)
        out.write('\t'.join(header_cols) + '\n')

    # read data, skipping invalid lines, and perform computations
    # that will result in new columns
    for i, line in enumerate(fh):
        total_lines += 1
        line = line.rstrip('\n')
        if not line or line.startswith('#'):
            skipped_lines += 1
            if not invalid_line:
                first_invalid_line = i + 1
                invalid_line = line
            continue
        fields = line.split('\t')
        typed_fields = None
        for fun, col_idx, mode, col_name, ex in ops:
            try:
                if typed_fields is None:
                    # Before performing the first column operation
                    # try to convert all fields of the current row from
                    # strings to their actual types (according to
                    # --column-types). If that attempt fails (e.g. because
                    # of a missing value in an integer column), simply
                    # copy the original fields unconverted and leave the
                    # situation to the outer exception handling according
                    # to the non-computable settings in effect.
                    try:
                        typed_fields = cast_types(*fields)
                    except Exception:
                        typed_fields = fields[:]
                        raise
                try:
                    new_val = fun(*typed_fields)
                except NameError as e:
                    # Python 3.10+ would have the problematic name
                    # available as e.name
                    if non_existent_col_pat.fullmatch(str(e)) and (
                        not args.fail_on_non_existent_columns
                    ):
                        # Looks like a reference to a non-existent column
                        # and we are not supposed to fail on it directly.
                        # Reraise and have it handled as a non-computable
                        # row.
                        raise
                    # NameErrors are not row-specific, but indicate a
                    # general problem with the user-supplied expression.
                    sys.exit(
                        'While parsing expression "%s" the following '
                        'problem occured: "%s"' % (ex, str(e))
                    )
            except Exception as e:
                if args.fail_on_non_computable:
                    sys.exit(
                        'Could not compute a new column value using "%s" on '
                        'line #%d: "%s".  Error was "%s"'
                        % (ex, i, line, str(e))
                    )
                if args.keep_non_computable:
                    # write the original line unchanged and stop computing
                    # for this line
                    out.write(line + '\n')
                    break
                if args.non_computable_blank:
                    new_val = ''
                elif args.non_computable_default is not None:
                    new_val = args.non_computable_default
                else:
                    # --skip_non_computable
                    # (which is default behavior, too)
                    # log that a line got skipped and why, then stop
                    # computing for this line
                    print(e)
                    skipped_lines += 1
                    if not invalid_line:
                        first_invalid_line = i + 1
                        invalid_line = line
                    break
            if mode is Mode.INSERT:
                fields.insert(col_idx, new_val)
                typed_fields.insert(col_idx, new_val)
            elif mode is Mode.REPLACE:
                fields[col_idx] = new_val
                typed_fields[col_idx] = new_val
            else:
                fields.append(new_val)
                typed_fields.append(new_val)
        else:
            fields = [
                format_float_positional(field)
                if args.avoid_scientific_notation and type(field) is float
                else str(field)
                for field in fields
            ]
            out.write('\t'.join(fields) + '\n')
            lines_computed += 1


valid_lines = total_lines - skipped_lines
if valid_lines > 0:
    print(
        'Computed new column values for %4.2f%% of %d lines written.'
        % (100.0 * lines_computed / valid_lines, valid_lines)
    )
elif args.fail_on_non_existent_columns:
    # Warn the user that there could be an issue with an expression.
    print(
        'Could not compute a new column for any input row!  '
        'Please check your expression(s) "%s" for problems.'
        % expr
    )
else:
    # Same, but the problem could also be a reference to a non-existent
    # column.
    print(
        'Could not compute a new column for any input row!  '
        'Please check your expression(s) "%s" for references to non-existent '
        'columns or other problems.'
        % expr
    )
if skipped_lines > 0:
    print('Skipped %d invalid lines starting at line #%d: "%s"' %
          (skipped_lines, first_invalid_line, invalid_line))
if lines_computed < valid_lines:
    print(
        'Rewrote %d lines unmodified because computation of a new value failed'
        % (valid_lines - lines_computed)
    )
