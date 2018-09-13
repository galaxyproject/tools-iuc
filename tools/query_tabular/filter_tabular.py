#!/usr/bin/env python

from __future__ import print_function

import json
import optparse
import os.path
import sys

from filters import filter_file


def __main__():
    # Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option('-i', '--input', dest='input', default=None,
                      help='Input file for filtering')
    parser.add_option('-j', '--jsonfile', dest='jsonfile', default=None,
                      help='JSON array of filter specifications')
    parser.add_option('-o', '--output', dest='output', default=None,
                      help='Output file for query results')
    parser.add_option('-v', '--verbose', dest='verbose', default=False,
                      action='store_true',
                      help='verbose')
    (options, args) = parser.parse_args()

    if options.input is not None:
        try:
            inputPath = os.path.abspath(options.input)
            inputFile = open(inputPath, 'r')
        except Exception as e:
            exit('Error: %s' % (e))
    else:
        inputFile = sys.stdin

    if options.output is not None:
        try:
            outputPath = os.path.abspath(options.output)
            outputFile = open(outputPath, 'w')
        except Exception as e:
            exit('Error: %s' % (e))
    else:
        outputFile = sys.stdout

    filters = None
    if options.jsonfile:
        try:
            with open(options.jsonfile) as fh:
                filters = json.load(fh)
        except Exception as e:
            exit('Error: %s' % (e))

    if options.verbose and filters:
        for f in filters:
            print('%s  %s' % (f['filter'],
                  ', '.join(
                  ['%s: %s' % (k, f[k])
                   for k in set(f.keys()) - set(['filter'])])),
                  file=sys.stdout)

    try:
        filter_file(inputFile, outputFile, filters=filters)
    except Exception as e:
        exit('Error: %s' % (e))


if __name__ == "__main__":
    __main__()
