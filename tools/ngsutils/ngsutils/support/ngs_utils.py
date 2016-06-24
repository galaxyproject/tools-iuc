#!/usr/bin/env python
"""
Common util classes / functions for the NGS project
"""
import collections
import gzip
import os
import re
import sys


def format_number(n):
    '''
    >>> format_number(1000)
    '1,000'
    >>> format_number(1234567)
    '1,234,567'
    '''
    ar = list(str(n))
    for i in range(len(ar))[::-3][1:]:
        ar.insert(i + 1, ',')
    return ''.join(ar)


def natural_sort(ar):
    '''
    >>> natural_sort('1 3 4 2 5'.split())
    ['1', '2', '3', '4', '5']
    >>> natural_sort('1 10 20 2 3 4'.split())
    ['1', '2', '3', '4', '10', '20']
    '''
    to_sort = []
    for item in ar:
        spl = re.split('(\d+)', item)
        l2 = []
        for el in spl:
            try:
                n = int(el)
            except:
                n = el
            l2.append(n)
        to_sort.append((l2, item))

    to_sort.sort()
    return [x[1] for x in to_sort]


def dictify(values, colnames):
    """
    Convert a list of values into a dictionary based upon given column names.

    If the column name starts with an '@', the value is assumed to be a comma
    separated list.

    If the name starts with a '#', the value is assumed to be an int.

    If the name starts with '@#', the value is assumed to  a comma separated
    list of ints.

    """
    d = {}
    for i in xrange(len(colnames)):
        key = colnames[i]
        split = False
        num = False

        if key[0] == '@':
            key = key[1:]
            split = True
        if key[0] == '#':
            key = key[1:]
            num = True

        if i < len(values):
            if num and split:
                val = [int(x) for x in values[i].rstrip(',').split(',')]
            elif num:
                val = int(values[i])
            elif split:
                val = values[i].rstrip(',').split(',')
            else:
                val = values[i]

            d[key] = val

        else:
            d[key] = None

    return d


def gzip_aware_open(fname):
    if fname == '-':
        f = sys.stdin
    elif fname[-3:] == '.gz' or fname[-4:] == '.bgz':
        f = gzip.open(os.path.expanduser(fname))
    else:
        f = open(os.path.expanduser(fname))
    return f


class gzip_opener:
    '''
    A Python 2.6 class to handle 'with' opening of text files that may
    or may not be gzip compressed.
    '''

    def __init__(self, fname):
        self.fname = fname

    def __enter__(self):
        self.f = gzip_aware_open(self.fname)
        return self.f

    def __exit__(self, type, value, traceback):
        if self.f != sys.stdin:
            self.f.close()
        return False


def filenames_to_uniq(names, new_delim='.'):
    '''
    Given a set of file names, produce a list of names consisting of the
    uniq parts of the names. This works from the end of the name.  Chunks of
    the name are split on '.' and '-'.

    For example:
        A.foo.bar.txt
        B.foo.bar.txt
        returns: ['A','B']

        AA.BB.foo.txt
        CC.foo.txt
        returns: ['AA.BB','CC']

    >>> filenames_to_uniq('a.foo.bar.txt b.foo.bar.txt'.split())
    ['a', 'b']
    >>> filenames_to_uniq('a.b.foo.txt c.foo.txt'.split())
    ['a.b', 'c']

    '''
    name_words = []
    maxlen = 0
    for name in names:
        name_words.append(name.replace('.', ' ').replace('-', ' ').strip().split())
        name_words[-1].reverse()
        if len(name_words[-1]) > maxlen:
            maxlen = len(name_words[-1])

    common = [False, ] * maxlen
    for i in xrange(maxlen):
        last = None
        same = True
        for nameword in name_words:
            if i >= len(nameword):
                same = False
                break
            if not last:
                last = nameword[i]
            elif nameword[i] != last:
                same = False
                break
        common[i] = same

    newnames = []
    for nameword in name_words:
        nn = []
        for (i, val) in enumerate(common):
            if not val and i < len(nameword):
                nn.append(nameword[i])
        nn.reverse()
        newnames.append(new_delim.join(nn))
    return newnames


def parse_args(argv, defaults=None, expected_argc=0):
    opts = {}
    if defaults:
        opts.update(defaults)

    args = []

    i = 0
    while i < len(argv):
        if argv[i][0] == '-':
            arg = argv[i].lstrip('-')
            if '=' in arg:
                k, v = arg.split('=', 2)
                if k in defaults:
                    if type(defaults[k]) == float:
                        opts[k] = float(v)
                    elif type(defaults[k]) == int:
                        opts[k] = int(v)
                    else:
                        opts[k] = v
            else:
                opts[arg] = True
        else:
            args.append(argv[i])
        i += 1

    while len(args) < expected_argc:
        args.append(None)
    return opts, args


class memoize(object):
    'Simple memoizing decorator to cache results'

    def __init__(self, func):
        self.func = func
        self.cache = {}

    def __call__(self, *args):
        if not isinstance(args, collections.Hashable):
            # uncacheable. a list, for instance.
            # better to not cache than blow up.
            return self.func(*args)

        if args in self.cache:
            return self.cache[args]
        else:
            value = self.func(*args)
            self.cache[args] = value
            return value
