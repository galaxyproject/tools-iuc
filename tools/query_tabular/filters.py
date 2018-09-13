#!/usr/binsenv python

from __future__ import print_function

import re
import sys


class LineFilter(object):
    def __init__(self, source, filter_dict):
        self.source = source
        self.filter_dict = filter_dict
        self.func = lambda i, l: l.rstrip('\r\n') if l else None
        self.src_lines = []
        self.src_line_cnt = 0
        if not filter_dict:
            return
        if filter_dict['filter'] == 'regex':
            rgx = re.compile(filter_dict['pattern'])
            if filter_dict['action'] == 'exclude_match':
                self.func = lambda i, l: l if not rgx.match(l) else None
            elif filter_dict['action'] == 'include_match':
                self.func = lambda i, l: l if rgx.match(l) else None
            elif filter_dict['action'] == 'exclude_find':
                self.func = lambda i, l: l if not rgx.search(l) else None
            elif filter_dict['action'] == 'include_find':
                self.func = lambda i, l: l if rgx.search(l) else None
        elif filter_dict['filter'] == 'select_columns':
            cols = [int(c) - 1 for c in filter_dict['columns']]
            self.func = lambda i, l: self.select_columns(l, cols)
        elif filter_dict['filter'] == 'replace':
            p = filter_dict['pattern']
            r = filter_dict['replace']
            c = int(filter_dict['column']) - 1
            self.func = lambda i, l: '\t'.join(
                [x if j != c else re.sub(p, r, x)
                 for j, x in enumerate(l.split('\t'))])
        elif filter_dict['filter'] == 'prepend_line_num':
            self.func = lambda i, l: '%d\t%s' % (i, l)
        elif filter_dict['filter'] == 'append_line_num':
            self.func = lambda i, l: '%s\t%d' % (l.rstrip('\r\n'), i)
        elif filter_dict['filter'] == 'prepend_text':
            s = filter_dict['column_text']
            self.func = lambda i, l: '%s\t%s' % (s, l)
        elif filter_dict['filter'] == 'append_text':
            s = filter_dict['column_text']
            self.func = lambda i, l: '%s\t%s' % (l.rstrip('\r\n'), s)
        elif filter_dict['filter'] == 'skip':
            cnt = filter_dict['count']
            self.func = lambda i, l: l if i > cnt else None
        elif filter_dict['filter'] == 'normalize':
            cols = [int(c) - 1 for c in filter_dict['columns']]
            sep = filter_dict['separator']
            self.func = lambda i, l: self.normalize(l, cols, sep)

    def __iter__(self):
        return self

    def __next__(self):
        if not self.src_lines:
            self.get_lines()
        if self.src_lines:
            return self.src_lines.pop(0)
        raise StopIteration

    next = __next__

    def select_columns(self, line, cols):
        fields = line.split('\t')
        return '\t'.join([fields[x] for x in cols])

    def normalize(self, line, split_cols, sep):
        lines = []
        fields = line.rstrip('\r\n').split('\t')
        split_fields = dict()
        cnt = 0
        for c in split_cols:
            if c < len(fields):
                split_fields[c] = fields[c].split(sep)
                cnt = max(cnt, len(split_fields[c]))
        if cnt == 0:
            lines.append('\t'.join(fields))
        else:
            for n in range(0, cnt):
                flds = [x if c not in split_cols else split_fields[c][n]
                        if n < len(split_fields[c])
                        else '' for (c, x) in enumerate(fields)]
                lines.append('\t'.join(flds))
        return lines

    def get_lines(self):
        for i, next_line in enumerate(self.source):
            self.src_line_cnt += 1
            line = self.func(self.src_line_cnt, next_line)
            if line:
                if isinstance(line, list):
                    self.src_lines.extend(line)
                else:
                    self.src_lines.append(line)
                return


class TabularReader:
    """
    Tabular file iterator. Returns a list
    """
    def __init__(self, input_file, skip=0, comment_char=None, col_idx=None,
                 filters=None):
        self.skip = skip
        self.comment_char = comment_char
        self.col_idx = col_idx
        self.filters = filters
        self.tsv_file = \
            input_file if hasattr(input_file, 'readline') else open(input_file)
        if skip and skip > 0:
            for i in range(skip):
                if not self.tsv_file.readline():
                    break
        source = LineFilter(self.tsv_file, None)
        if comment_char:
            source = LineFilter(source,
                                {"filter": "regex", "pattern": comment_char,
                                 "action": "exclude_match"})
        if filters:
            for f in filters:
                source = LineFilter(source, f)
        self.source = source

    def __iter__(self):
        return self

    def __next__(self):
        ''' Iteration '''
        for i, line in enumerate(self.source):
            fields = line.rstrip('\r\n').split('\t')
            if self.col_idx:
                fields = [fields[i] for i in self.col_idx]
            return fields
        raise StopIteration

    next = __next__


def filter_file(input_file, output, skip=0, comment_char='#', filters=None):
    data_lines = 0
    try:
        tr = TabularReader(input_file, skip=skip, comment_char=comment_char,
                           filters=filters)
        for linenum, fields in enumerate(tr):
            data_lines += 1
            try:
                output.write('%s\n' % '\t'.join(fields))
            except Exception as e:
                print('Failed at line: %d err: %s' % (linenum, e),
                      file=sys.stderr)
    except Exception as e:
        exit('Error: %s' % (e))
