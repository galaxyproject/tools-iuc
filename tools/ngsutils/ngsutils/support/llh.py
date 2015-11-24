'''
Methods for calculating log-likelihoods for nucleotide frequencies
'''
import math
import collections
from ngsutils.support import memoize

_default_background = {'A': 0.3, 'T': 0.3, 'C': 0.2, 'G': 0.2}

NucleotideLogLikelihood = collections.namedtuple('NucleotideLogLikelihood', 'A C G T pseudo')

@memoize
def pseudo_count(N, bg):
    '''
    >>> pseudo_count(100, _default_background['A'])
    3
    >>> pseudo_count(100, _default_background['C'])
    2
    '''

    return bg * math.sqrt(N)


def calc_llh(A, C, G, T, bg=_default_background, pseudo='auto'):
    if pseudo == 'auto':
        N = A + C + G + T
        Ap = float(A) + pseudo_count(N, bg['A'])
        Cp = float(C) + pseudo_count(N, bg['C'])
        Gp = float(G) + pseudo_count(N, bg['G'])
        Tp = float(T) + pseudo_count(N, bg['T'])
    elif pseudo:
        Ap = float(A) + pseudo
        Cp = float(C) + pseudo
        Gp = float(G) + pseudo
        Tp = float(T) + pseudo
    else:
        Ap = float(A)
        Cp = float(C)
        Gp = float(G)
        Tp = float(T)

    Np = Ap + Cp + Gp + Tp

    freqA = Ap / Np
    freqC = Cp / Np
    freqG = Gp / Np
    freqT = Tp / Np

    return NucleotideLogLikelihood(math.log(freqA / bg['A']), math.log(freqC / bg['C']), math.log(freqG / bg['G']), math.log(freqT / bg['T']), pseudo)



if __name__ == '__main__':
    import doctest
    doctest.testmod()
