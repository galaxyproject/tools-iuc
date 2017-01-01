import collections
import gzip
import os
import re
import sys
try:
    from eta import ETA
except:
    pass


class FASTARead(collections.namedtuple('FASTARecord', 'name comment seq')):
    def __repr__(self):
        if self.comment:
            return '>%s %s\n%s\n' % (self.name, self.comment, self.seq)
        return '>%s\n%s\n' % (self.name, self.seq)

    def subseq(self, start, end, comment=None):
        if self.comment:
            comment = '%s %s' % (self.comment, comment)

        return FASTARead(self.name, comment, self.seq[start:end])

    def clone(self, name=None, comment=None, seq=None):
        n = name if name else self.name
        c = comment if comment else self.comment
        s = seq if seq else self.seq

        return FASTARead(n, c, s)

    def write(self, out):
        out.write(repr(self))


class FASTA(object):
    def __init__(self, fname=None, fileobj=None, qual=False):
        self.fname = fname
        self.qual = qual
        if fileobj:
            self.fileobj = fileobj
        else:
            if self.fname == '-':
                self.fileobj = sys.stdin
            elif self.fname[-3:] == '.gz' or self.fname[-4:] == '.bgz':
                self.fileobj = gzip.open(os.path.expanduser(self.fname))
            else:
                self.fileobj = open(os.path.expanduser(self.fname))

        if not self.fileobj:
            raise ValueError("Missing valid filename or fileobj")

    def close(self):
        if self.fileobj != sys.stdout:
            self.fileobj.close()

    def tell(self):
        # always relative to uncompressed...
        return self.fileobj.tell()

    def seek(self, pos, whence=0):
        self.fileobj.seek(pos, whence)

    def fetch(self, quiet=False):
        name = ''
        comment = ''
        seq = ''

        if not quiet and self.fname and self.fname != '-':
            eta = ETA(os.stat(self.fname).st_size, fileobj=self.fileobj)
        else:
            eta = None

        for line in self.fileobj:
            line = line.strip()
            if not line:
                continue
            if line[0] == '#':
                continue

            if line[0] == '>':
                if name and seq:
                    if eta:
                        eta.print_status(extra=name)
                    yield FASTARead(name, comment, seq)

                spl = re.split(r'[ \t]', line[1:], maxsplit=1)
                name = spl[0]
                if len(spl) > 1:
                    comment = spl[1]
                else:
                    comment = ''
                seq = ''

            else:
                if self.qual:
                    seq = seq + ' ' + line
                else:
                    seq += line

        if name and seq:
            if eta:
                eta.print_status(extra=name)
            yield FASTARead(name, comment, seq)

        if eta:
            eta.done()


def gzip_reader(fname, quiet=False, callback=None, done_callback=None, fileobj=None):
    if fileobj:
        f = fileobj
    elif fname == '-':
        f = sys.stdin
    elif fname[-3:] == '.gz' or fname[-4:] == '.bgz':
        f = gzip.open(os.path.expanduser(fname))
    else:
        f = open(os.path.expanduser(fname))

    if quiet or fname == '-':
        eta = None
    else:
        eta = ETA(os.stat(fname).st_size, fileobj=f)

    for line in f:
        if eta:
            if callback:
                extra = callback()
            else:
                extra = ''

            eta.print_status(extra=extra)
        yield line

        if done_callback and done_callback():
                break

    if f != sys.stdin:
        f.close()

    if eta:
        eta.done()


class Symbolize(object):
    'Converts strings to symbols - basically a cache of strings'

    def __init__(self):
        self.__cache = {}

    def __getitem__(self, k):
        if k not in self.__cache:
            self.__cache[k] = k

        return self.__cache[k]


symbols = Symbolize()

_compliments = {
    'a': 't',
    'A': 'T',
    'c': 'g',
    'C': 'G',
    'g': 'c',
    'G': 'C',
    't': 'a',
    'T': 'A',
    'n': 'n',
    'N': 'N'
}


def revcomp(seq):
    '''
    >>> revcomp('ATCGatcg')
    'cgatCGAT'
    '''
    ret = []

    for s in seq:
        ret.append(_compliments[s])

    ret.reverse()
    return ''.join(ret)


class Counts(object):
    '''
    Setup simple binning.  Bins are continuous 0->max.  Values are added to
    bins and then means / distributions can be calculated.
    '''

    def __init__(self):
        self.bins = []

    def add(self, val):
        while len(self.bins) <= val:
            self.bins.append(0)
        self.bins[val] += 1

    def mean(self):
        acc = 0
        count = 0

        for i, val in enumerate(self.bins):
            acc += (i * val)
            count += val

        if count > 0:
            return float(acc) / count

    def max(self):
        return len(self.bins) - 1


def memoize(func):
    if 'TESTING' in os.environ or 'DEBUG' in os.environ:
        return func

    __cache = {}

    def inner(*args, **kwargs):
        k = (args, tuple(kwargs.iteritems()))
        if k not in __cache:
            __cache[k] = func(*args, **kwargs)
        return __cache[k]

    inner.__doc__ = '(@memoized %s)\n%s' % (func.__name__, func.__doc__)
    return inner


def quoted_split(s, delim, quote_char='"'):
    tokens = []

    buf = ""
    inquote = False

    for c in s:
        if inquote:
            buf += c
            if c == quote_char:
                inquote = False
        elif c == delim:
            tokens.append(buf)
            buf = ""
        else:
            buf += c
            if c == quote_char:
                inquote = True

    if buf:
        tokens.append(buf)

    return tokens
