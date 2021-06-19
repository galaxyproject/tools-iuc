#!/usr/bin/env python

from __future__ import print_function

import math
import re
import sqlite3 as sqlite
import sys


TABLE_QUERY = \
    """
    SELECT name, sql
    FROM sqlite_master
    WHERE type='table'
    ORDER BY name
    """


def msg(e):
    print(e, file=sys.stderr)


def regex_match(expr, item):
    return re.match(expr, item) is not None


def regex_search(expr, item):
    return re.search(expr, item) is not None


def regex_sub(expr, replace, item):
    return re.sub(expr, replace, item)


def math_acos(x):
    try:
        return math.acos(x)
    except ValueError as ve:
        msg('acos(%s): %s' % (x, ve))
        return None


def math_acosh(x):
    try:
        return math.acosh(x)
    except ValueError as ve:
        msg(f'acosh({x}): {ve}')
        return None


def math_asin(x):
    try:
        return math.asin(x)
    except ValueError as ve:
        msg(f'asin({x}): {ve}')
        return None


def math_asinh(x):
    try:
        return math.asinh(x)
    except ValueError as ve:
        msg(f'asinh({x}): {ve}')
        return None


def math_atan(x):
    try:
        return math.atan(x)
    except ValueError as ve:
        msg(f'atan({x}): {ve}')
        return None


def math_atanh(x):
    try:
        return math.atanh(x)
    except ValueError as ve:
        msg(f'atanh({x}): {ve}')
        return None


def math_atan2(x, y):
    try:
        return math.atan2(x, y)
    except ValueError as ve:
        msg(f'atan2({x}, {y}): {ve}')
        return None


def math_ceil(x):
    try:
        return math.ceil(x)
    except ValueError as ve:
        msg(f'ceil({x}): {ve}')
        return None


def math_cos(x):
    try:
        return math.cos(x)
    except ValueError as ve:
        msg(f'cos({x}): {ve}')
        return None


def math_cosh(x):
    try:
        return math.cosh(x)
    except ValueError as ve:
        msg(f'cosh({x}): {ve}')
        return None


def math_degrees(x):
    try:
        return math.degrees(x)
    except ValueError as ve:
        msg(f'degrees({x}): {ve}')
        return None


def math_exp(x):
    try:
        return math.exp(x)
    except ValueError as ve:
        msg(f'exp({x}): {ve}')
        return None


def math_expm1(x):
    try:
        return math.expm1(x)
    except ValueError as ve:
        msg(f'expm1({x}): {ve}')
        return None


def math_fabs(x):
    try:
        return math.fabs(x)
    except ValueError as ve:
        msg(f'fabs({x}): {ve}')
        return None


def math_floor(x):
    try:
        return math.floor(x)
    except ValueError as ve:
        msg(f'floor({x}): {ve}')
        return None


def math_fmod(x, y):
    try:
        return math.fmod(x, y)
    except ValueError as ve:
        msg(f'fmod({x}, {y}): {ve}')
        return None


def math_blog(b, x):
    try:
        return math.log(b, x)
    except ValueError as ve:
        msg(f'log({b}, {x}): {ve}')
        return None


def math_log(x):
    try:
        return math.log(x)
    except ValueError as ve:
        msg(f'log({x}): {ve}')
        return None


def math_log10(x):
    try:
        return math.log10(x)
    except ValueError as ve:
        msg(f'log10({x}): {ve}')
        return None


def math_log1p(x):
    try:
        return math.log1p(x)
    except ValueError as ve:
        msg(f'log1p({x}): {ve}')
        return None


def math_log2(x):
    try:
        return math.log2(x)
    except ValueError as ve:
        msg(f'log2({x}): {ve}')
        return None


def math_mod(x, y):
    try:
        return x % y
    except ValueError as ve:
        msg(f'mod({x}, {y}): {ve}')
        return None


def math_pow(x, y):
    try:
        return math.pow(x, y)
    except ValueError as ve:
        msg(f'pow({x}, {y}): {ve}')
        return None


def math_radians(x):
    try:
        return math.radians(x)
    except ValueError as ve:
        msg(f'radians({x}): {ve}')
        return None


def math_sin(x):
    try:
        return math.sin(x)
    except ValueError as ve:
        msg(f'sin({x}): {ve}')
        return None


def math_sinh(x):
    try:
        return math.sinh(x)
    except ValueError as ve:
        msg(f'sinh({x}): {ve}')
        return None


def math_sqrt(x):
    try:
        return math.sqrt(x)
    except ValueError as ve:
        msg(f'sqrt({x}): {ve}')
        return None


def math_tan(x):
    try:
        return math.tan(x)
    except ValueError as ve:
        msg(f'tan({x}): {ve}')
        return None


def math_tanh(x):
    try:
        return math.tanh(x)
    except ValueError as ve:
        msg(f'tanh({x}): {ve}')
        return None


def math_trunc(x):
    try:
        return math.trunc(x)
    except ValueError as ve:
        msg(f'trunc({x}): {ve}')
        return None


def get_connection(sqlitedb_path, addfunctions=True):
    sqlite.enable_callback_tracebacks(addfunctions)
    conn = sqlite.connect(sqlitedb_path)
    if addfunctions:
        conn.create_function("re_match", 2, regex_match)
        conn.create_function("re_search", 2, regex_search)
        conn.create_function("re_sub", 3, regex_sub)
        conn.create_function("acos", 1, math_acos)
        conn.create_function("acosh", 1, math_acosh)
        conn.create_function("asin", 1, math_asin)
        conn.create_function("asinh", 1, math_asinh)
        conn.create_function("atan", 1, math_atan)
        conn.create_function("atanh", 1, math_atanh)
        conn.create_function("atan2", 2, math_atan2)
        conn.create_function("ceil", 1, math_ceil)
        conn.create_function("cos", 1, math_cos)
        conn.create_function("cosh", 1, math_cosh)
        conn.create_function("degrees", 1, math_degrees)
        conn.create_function("exp", 1, math_exp)
        conn.create_function("expm1", 1, math_expm1)
        conn.create_function("fabs", 1, math_fabs)
        conn.create_function("floor", 1, math_floor)
        conn.create_function("fmod", 2, math_fmod)
        conn.create_function("log", 1, math_log)
        conn.create_function("log", 2, math_blog)
        conn.create_function("log10", 1, math_log10)
        conn.create_function("log2", 1, math_log2)
        conn.create_function("log1p", 1, math_log1p)
        conn.create_function("mod", 2, math_mod)
        conn.create_function("pow", 2, math_pow)
        conn.create_function("radians", 1, math_radians)
        conn.create_function("sin", 1, math_sin)
        conn.create_function("sinh", 1, math_sinh)
        conn.create_function("sqrt", 1, math_sqrt)
        conn.create_function("tan", 1, math_tan)
        conn.create_function("tanh", 1, math_tanh)
        conn.create_function("trunc", 1, math_trunc)
    return conn


def describe_tables(conn, outputFile):
    try:
        c = conn.cursor()
        tables_query = TABLE_QUERY
        rslt = c.execute(tables_query).fetchall()
        for table, sql in rslt:
            print("Table %s:" % table, file=outputFile)
            try:
                col_query = 'SELECT * FROM %s LIMIT 0' % table
                cur = conn.cursor().execute(col_query)
                cols = [col[0] for col in cur.description]
                print(" Columns: %s" % cols, file=outputFile)
            except Exception as exc:
                print("Warning: %s" % exc, file=sys.stderr)
    except Exception as e:
        exit('describe_tables Error: %s' % (e))
    exit(0)


def run_query(conn, query, outputFile, no_header=False, comment_char='#'):
    cur = conn.cursor()
    results = cur.execute(query)
    if outputFile is not None:
        if not no_header:
            outputFile.write("%s%s\n" % (comment_char, '\t'.join(
                [str(col[0]) for col in cur.description])))
        for i, row in enumerate(results):
            outputFile.write("%s\n" % '\t'.join(
                [str(val) if val is not None else '' for val in row]))
    else:
        conn.commit()
        return results
