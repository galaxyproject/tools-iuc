#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2011-2015, Michigan State University.
# Copyright (C) 2015, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
from __future__ import print_function

import os
import sys

import screed
from khmer import Countgraph, ReadParser
from khmer.utils import (broken_paired_reader, write_record)

CUTOFF = 50


def main():
    counting_ht = sys.argv[1]
    infiles = sys.argv[2:]

    print('file with ht: %s' % counting_ht)

    print('making hashtable')
    ht = Countgraph.load(counting_ht)
    K = ht.ksize()

    for infile in infiles:
        print('filtering', infile)
        outfile = os.path.basename(infile) + '.below'

        outfp = open(outfile, 'w')

        paired_iter = broken_paired_reader(ReadParser(infile), min_length=K,
                                           force_single=True)
        for n, is_pair, read1, read2 in paired_iter:
            name = read1.name
            seq = read1.sequence
            if 'N' in seq:
                return None, None

            trim_seq, trim_at = ht.trim_below_abundance(seq, CUTOFF)

            if trim_at >= K:
                write_record(screed.Record(name=name, sequence=trim_seq), outfp)


if __name__ == '__main__':
    main()
