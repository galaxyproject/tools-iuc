#!/usr/bin/env python

from __future__ import print_function

import re
import sys

from filters import TabularReader


SQLITE_KEYWORDS = [
    'ABORT',
    'ACTION',
    'ADD',
    'AFTER',
    'ALL',
    'ALTER',
    'ANALYZE',
    'AND',
    'AS',
    'ASC',
    'ATTACH',
    'AUTOINCREMENT',
    'BEFORE',
    'BEGIN',
    'BETWEEN',
    'BY',
    'CASCADE',
    'CASE',
    'CAST',
    'CHECK',
    'COLLATE',
    'COLUMN',
    'COMMIT',
    'CONFLICT',
    'CONSTRAINT',
    'CREATE',
    'CROSS',
    'CURRENT_DATE',
    'CURRENT_TIME',
    'CURRENT_TIMESTAMP',
    'DATABASE',
    'DEFAULT',
    'DEFERRABLE',
    'DEFERRED',
    'DELETE',
    'DESC',
    'DETACH',
    'DISTINCT',
    'DROP',
    'EACH',
    'ELSE',
    'END',
    'ESCAPE',
    'EXCEPT',
    'EXCLUSIVE',
    'EXISTS',
    'EXPLAIN',
    'FAIL',
    'FOR',
    'FOREIGN',
    'FROM',
    'FULL',
    'GLOB',
    'GROUP',
    'HAVING',
    'IF',
    'IGNORE',
    'IMMEDIATE',
    'IN',
    'INDEX',
    'INDEXED',
    'INITIALLY',
    'INNER',
    'INSERT',
    'INSTEAD',
    'INTERSECT',
    'INTO',
    'IS',
    'ISNULL',
    'JOIN',
    'KEY',
    'LEFT',
    'LIKE',
    'LIMIT',
    'MATCH',
    'NATURAL',
    'NO',
    'NOT',
    'NOTNULL',
    'NULL',
    'OF',
    'OFFSET',
    'ON',
    'OR',
    'ORDER',
    'OUTER',
    'PLAN',
    'PRAGMA',
    'PRIMARY',
    'QUERY',
    'RAISE',
    'RECURSIVE',
    'REFERENCES',
    'REGEXP',
    'REINDEX',
    'RELEASE',
    'RENAME',
    'REPLACE',
    'RESTRICT',
    'RIGHT',
    'ROLLBACK',
    'ROW',
    'SAVEPOINT',
    'SELECT',
    'SET',
    'TABLE',
    'TEMP',
    'TEMPORARY',
    'THEN',
    'TO',
    'TRANSACTION',
    'TRIGGER',
    'UNION',
    'UNIQUE',
    'UPDATE',
    'USING',
    'VACUUM',
    'VALUES',
    'VIEW',
    'VIRTUAL',
    'WHEN',
    'WHERE',
    'WITH',
    'WITHOUT'
]


def get_valid_column_name(name):
    valid_name = name
    if not name or not name.strip():
        return None
    elif name.upper() in SQLITE_KEYWORDS:
        valid_name = '"%s"' % name
    elif re.match('^[a-zA-Z]\w*$', name):
        pass
    elif re.match('^"[^"]+"$', name):
        pass
    elif re.match('^\[[^\[\]]*\]$', name):
        pass
    elif re.match("^`[^`]+`$", name):
        pass
    elif name.find('"') < 0:
        valid_name = '"%s"' % name
    elif name.find('[') < 0 and name.find(']') < 0:
        valid_name = '[%s]' % name
    elif name.find('`') < 0:
        valid_name = '`%s`' % name
    elif name.find("'") < 0:
        valid_name = "'%s'" % name
    return valid_name


def getValueType(val):
    if val or 0. == val:
        try:
            int(val)
            return 'INTEGER'
        except ValueError:
            try:
                float(val)
                return 'REAL'
            except ValueError:
                return 'TEXT'
    return None


def get_column_def(file_path, table_name, skip=0, comment_char='#',
                   column_names=None, max_lines=100, load_named_columns=False,
                   firstlinenames=False, filters=None):
    col_pref = ['TEXT', 'REAL', 'INTEGER', None]
    col_types = []
    col_idx = None
    col_names = []
    try:
        tr = TabularReader(file_path, skip=skip, comment_char=comment_char,
                           col_idx=None, filters=filters)
        for linenum, fields in enumerate(tr):
            if linenum == 0 and firstlinenames:
                col_names = [get_valid_column_name(name) or 'c%d' % (i + 1)
                             for i, name in enumerate(fields)]
                continue
            if linenum > max_lines:
                break
            try:
                while len(col_types) < len(fields):
                    col_types.append(None)
                for i, val in enumerate(fields):
                    colType = getValueType(val)
                    if col_pref.index(colType) < col_pref.index(col_types[i]):
                        col_types[i] = colType
            except Exception as e:
                print('Failed at line: %d err: %s' % (linenum, e),
                      file=sys.stderr)
    except Exception as e:
        print('Failed: %s' % (e), file=sys.stderr)
    for i, col_type in enumerate(col_types):
        if not col_type:
            col_types[i] = 'TEXT'
    if not col_names:
        col_names = ['c%d' % i for i in range(1, len(col_types) + 1)]
    if column_names:
        if load_named_columns:
            col_idx = []
            cnames = []
            for i, cname in enumerate(
                    [cn.strip() for cn in column_names.split(',')]):
                if cname != '':
                    col_idx.append(i)
                    cnames.append(cname)
            col_types = [col_types[i] for i in col_idx]
            col_names = cnames
        else:
            for i, cname in enumerate(
                    [cn.strip() for cn in column_names.split(',')]):
                if cname and i < len(col_names):
                    col_names[i] = cname
    col_def = []
    for i, col_name in enumerate(col_names):
        col_def.append('%s %s' % (col_names[i], col_types[i]))
    return col_names, col_types, col_def, col_idx


def create_table(conn, file_path, table_name, skip=0, comment_char='#',
                 pkey_autoincr=None, column_names=None,
                 load_named_columns=False, firstlinenames=False,
                 filters=None, unique_indexes=[], indexes=[]):
    col_names, col_types, col_def, col_idx = \
        get_column_def(file_path, table_name, skip=skip,
                       comment_char=comment_char, column_names=column_names,
                       load_named_columns=load_named_columns,
                       firstlinenames=firstlinenames,
                       filters=filters)
    col_func = [float if t == 'REAL' else int
                if t == 'INTEGER' else str for t in col_types]
    table_def = 'CREATE TABLE %s (\n    %s%s\n);' % (
                table_name,
                '%s INTEGER PRIMARY KEY AUTOINCREMENT,' %
                pkey_autoincr if pkey_autoincr else '',
                ', \n    '.join(col_def))
    # print >> sys.stdout, table_def
    insert_stmt = 'INSERT INTO %s(%s) VALUES(%s)' % (
                  table_name, ','.join(col_names),
                  ','.join(["?" for x in col_names]))
    # print >> sys.stdout, insert_stmt
    data_lines = 0
    try:
        c = conn.cursor()
        c.execute(table_def)
        conn.commit()
        c.close()

        c = conn.cursor()
        tr = TabularReader(file_path, skip=skip, comment_char=comment_char,
                           col_idx=col_idx, filters=filters)
        for linenum, fields in enumerate(tr):
            if linenum == 0 and firstlinenames:
                continue
            data_lines += 1
            try:
                vals = [col_func[i](x)
                        if x else None for i, x in enumerate(fields)]
                c.execute(insert_stmt, vals)
            except Exception as e:
                print('Failed at line: %d err: %s' % (linenum, e),
                      file=sys.stderr)
        conn.commit()
        c.close()
        for i, index in enumerate(unique_indexes):
            index_name = 'idx_uniq_%s_%d' % (table_name, i)
            index_columns = index.split(',')
            create_index(conn, table_name, index_name, index_columns,
                         unique=True)
        for i, index in enumerate(indexes):
            index_name = 'idx_%s_%d' % (table_name, i)
            index_columns = index.split(',')
            create_index(conn, table_name, index_name, index_columns)
    except Exception as e:
        exit('Error: %s' % (e))


def create_index(conn, table_name, index_name, index_columns, unique=False):
    index_def = 'CREATE %s INDEX %s on %s(%s)' % (
                'UNIQUE' if unique else '', index_name,
                table_name, ','.join(index_columns))
    try:
        c = conn.cursor()
        c.execute(index_def)
        conn.commit()
        c.close()
    except Exception as e:
        print('Failed: %s err: %s' % (index_def, e), file=sys.stderr)
        raise(e)
