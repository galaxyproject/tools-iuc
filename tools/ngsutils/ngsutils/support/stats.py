'''
various statistical tests and methods...
'''
import math

from ngsutils.support import memoize


def median(vals):
    '''
    >>> median([1,2,3])
    2
    >>> median([1,2,3,4])
    2.5
    '''
    vals.sort()

    if len(vals) % 2 == 1:
        return vals[len(vals) / 2]
    else:
        a = vals[(len(vals) / 2) - 1]
        b = vals[(len(vals) / 2)]
        return float(a + b) / 2


def mean_stdev(l):
    '''
    >>> mean_stdev([1,2,2,2])
    (1.75, 0.5)
    >>> mean_stdev([2,2,2,2])
    (2.0, 0.0)
    '''

    acc = 0
    for el in l:
        acc += el

    mean = float(acc) / len(l)
    acc = 0
    for el in l:
        acc += (el - mean) ** 2

    if len(l) > 2:
        stdev = math.sqrt(float(acc) / (len(l) - 1))
    else:
        stdev = 0.0

    return (mean, stdev)


def counts_median(d):
    '''
    Calculate the median from counts stored in a dictionary
    >>> counts_median({ 1: 4, 2: 1, 3: 4 })
    2
    >>> counts_median({ 1: 4, 3: 4 })
    2

    '''
    count = 0
    for k in d:
        count += d[k]

    if count == 0:
        return 0

    acc = 0.0
    last = 0
    for k in sorted(d):
        if last:
            return (last + k) / 2
        acc += d[k]
        if acc / count == 0.5:
            last = k
        elif acc / count > 0.5:
            return k


def counts_mean_stdev(d):
    '''

    calc mean / stdev when data is stored as counts in a dictionary

    Ex:
        { 1: 4, 2: 1, 3: 4 } = [1, 1, 1, 1, 2, 3, 3, 3, 3]

    >>> counts_mean_stdev({ 1: 4, 2: 1, 3: 4 })
    (2.0, 1.0)

    '''

    acc = 0
    count = 0
    for k in d:
        acc += k * d[k]
        count += d[k]

    mean = float(acc) / count

    acc = 0
    for k in d:
        acc += (((k - mean) ** 2) * d[k])

    if count > 2:
        stdev = math.sqrt(float(acc) / (count - 1))
    else:
        stdev = 0.0

    return (mean, stdev)


@memoize
def poisson_prob(x, mean):
    '''
        Return the probability that you could get x counts in
        a Poisson test with a mean value.

        prob(x) = sum(i=1..x){poisson(i)}

        >>> poisson_prob(6,10)
        0.1300960209527205
        >>> poisson_prob(8,10)
        0.33277427882095645
    '''
    acc = 0.0
    for i in xrange(1, x + 1):
        acc += poisson_func(i, mean)
    return acc


@memoize
def poisson_func(mu, lambd):
    '''
        This is the Poisson distribution function

        p(mu) = (lambda^mu * e^(-lambda)) / (mu!)

        mu is a count
        lambd is the mean

        >>> poisson_func(1,10)
        0.00045399929762484856
        >>> poisson_func(2,10)
        0.0022699964881242427
        >>> poisson_func(3,10)
        0.007566654960414142
    '''
    return (lambd ** mu) * (math.exp(-1 * lambd)) / _factorial(mu)


@memoize
def _factorial(x):
    '''
    >>> _factorial(1)
    1
    >>> _factorial(2)
    2
    >>> _factorial(3)
    6
    '''
    return math.factorial(x)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
