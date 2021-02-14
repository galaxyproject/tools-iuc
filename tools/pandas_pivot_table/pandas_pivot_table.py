#!/usr/bin/env python

import argparse
import json
import sys
from json.decoder import JSONDecodeError

import pandas as pd


def __main__():
    p = argparse.ArgumentParser()
    p.add_argument(
        '-i', '--input',
        type=argparse.FileType('r'),
        required=True,
        help='Tabular input file to pivot'
    )
    p.add_argument(
        '-o', '--output',
        type=argparse.FileType('w'),
        required=True,
        help='Output file'
    )
    p.add_argument(
        '-S', '--skiprows',
        type=int,
        default=0,
        help='Input column names'
    )
    p.add_argument(
        '-H', '--header',
        default=None,
        help='Input column names'
    )
    p.add_argument(
        '-P', '--prefix',
        default=None,
        help='Prefix for input column names'
    )
    p.add_argument(
        '-I', '--index',
        help='index columns'
    )
    p.add_argument(
        '-C', '--columns',
        help='columns values which are returned as columns'
    )
    p.add_argument(
        '-V', '--values',
        help='values'
    )
    p.add_argument(
        '-F', '--aggfunc',
        help='aggregate functions on the values'
    )
    p.add_argument(
        '-N', '--fill_value',
        default=None,
        help='fill value for missing values'
    )
    p.add_argument(
        '-f', '--float_format',
        default='%0.6f',
        help=''
    )
    args = p.parse_args()

    def getValueType(val):
        if val or 0. == val:
            try:
                return int(val)
            except ValueError:
                try:
                    return float(val)
                except ValueError:
                    return val
        return None

    def getColumn(name, dfcols, value_cols=None):
        dfname = None
        if name in dfcols:
            dfname = name
        else:
            try:
                i = int(name) - 1
                dfname = dfcols[i]
            except IndexError:
                sys.exit('%s not an index into %s' % (name, dfcols))
            except ValueError:
                sys.exit('%s not a column in %s' % (name, dfcols))
        if value_cols and dfname not in value_cols:
            sys.exit('%s not a value column in %s' % (name, value_cols))
        return dfname

    def getColumns(val, dfcols):
        fields = [v.strip() for v in val.split(',')]
        cols = []
        for name in fields:
            cols.append(getColumn(name, dfcols))
        return cols

    def getAggFunc(funcStr, dfcols, value_cols):
        af = funcStr
        try:
            af = json.loads(funcStr)
        except JSONDecodeError as de:
            sys.exit('"%s" is not a json string: %s' % (funcStr, de.msg))
        if isinstance(af, dict):
            aggfunc = {getColumn(k, dfcols, value_cols): v
                       for k, v in af.items()}
        elif isinstance(af, list):
            aggfunc = af
        else:
            aggfunc = af
        return aggfunc

    if args.prefix:
        df = pd.read_table(args.input,
                           skiprows=args.skiprows,
                           header=None,
                           prefix=args.prefix)
    elif args.header:
        df = pd.read_table(args.input,
                           skiprows=args.skiprows,
                           header=args.header)
    else:
        df = pd.read_table(args.input, skiprows=args.skiprows)
    df_columns = df.columns.tolist()
    index = getColumns(args.index, df_columns)
    columns = getColumns(args.columns, df_columns)
    values = getColumns(args.values, df_columns)
    fill_value = getValueType(args.fill_value)
    aggfunc = getAggFunc(args.aggfunc.replace('\'', '"'), df_columns, values)
    pdf = df.pivot_table(index=index, columns=columns,
                         values=values, aggfunc=aggfunc,
                         fill_value=fill_value)
    pdf_cols = ['_'.join([str(x) for x in reversed(p)])
                if isinstance(p, tuple) else str(p)
                for p in pdf.columns.tolist()]
    pdf.to_csv(args.output,
               sep='\t',
               float_format=args.float_format,
               header=pdf_cols)


if __name__ == "__main__":
    __main__()
