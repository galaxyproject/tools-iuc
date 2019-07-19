#!/usr/bin/env python

# Reports a beta diversity matrix for tabular input file
# using scikit-bio
# Daniel Blankenberg
from __future__ import print_function
import codecs
import optparse
import sys

from skbio import TreeNode
from skbio.diversity import beta_diversity


__VERSION__ = "0.0.1"

DELIMITER = '\t'

NEEDS_TREE = [ 'unweighted_unifrac', 'weighted_unifrac' ]

NEEDS_OTU_NAMES = [ 'unweighted_unifrac', 'weighted_unifrac' ]


def __main__():
    parser = optparse.OptionParser( usage="%prog [options]" )
    parser.add_option( '-v', '--version', dest='version', action='store_true', default=False, help='print version and exit' )
    parser.add_option( '-i', '--input', dest='input', action='store', type="string", default=None, help='Input abundance Filename' )
    parser.add_option( '', '--otu_column', dest='otu_column', action='store', type="int", default=None, help='OTU ID Column (1 based)' )
    parser.add_option( '', '--sample_columns', dest='sample_columns', action='store', type="string", default=None, help='Comma separated list of sample columns, unset to use all.' )
    parser.add_option( '', '--header', dest='header', action='store_true', default=False, help='Abundance file has a header line' )
    parser.add_option( '', '--distance_metric', dest='distance_metric', action='store', type="string", default=None, help='Distance metric to use' )
    parser.add_option( '', '--tree', dest='tree', action='store', type="string", default=None, help='Newick Tree Filename' )
    parser.add_option( '-o', '--output', dest='output', action='store', type="string", default=None, help='Output Filename' )
    (options, args) = parser.parse_args()
    if options.version:
        print("scikit-bio betadiversity from tabular file", __VERSION__, file=sys.stderr)
        sys.exit()

    if options.otu_column is not None:
        otu_column = options.otu_column - 1
    else:
        otu_column = None

    if options.sample_columns is None:
        with open( options.input, 'rb' ) as fh:
            line = fh.readline()
            columns = list(range( len( line.split( DELIMITER ) )))
            if otu_column in columns:
                columns.remove( otu_column )
    else:
        columns = [int( x ) - 1 for x in options.sample_columns.split( "," )]

    max_col = max( columns + [otu_column] )
    counts = [ [] for x in columns ]
    sample_names = []
    otu_names = []
    with open( options.input, 'rb' ) as fh:
        if options.header:
            header = fh.readline().rstrip('\n\r').split( DELIMITER )
            sample_names = [ header[i] for i in columns ]
        else:
            sample_names = [ "SAMPLE_%i" % x for x in range( len( columns ) ) ]
        for i, line in enumerate( fh ):
            fields = line.rstrip('\n\r').split( DELIMITER )
            if len(fields) <= max_col:
                print("Bad data line: ", fields, file=sys.stederr)
                continue
            if otu_column is not None:
                otu_names.append( fields[ otu_column ] )
            else:
                otu_names.append( "OTU_%i" % i )
            for j, col in enumerate( columns ):
                counts[ j ].append( int( fields[ col ] ) )

    extra_kwds = {}
    if options.distance_metric in NEEDS_OTU_NAMES:
        extra_kwds['otu_ids'] = otu_names
    if options.distance_metric in NEEDS_TREE:
        assert options.tree, Exception( "You must provide a newick tree when using '%s'" % options.distance_metric )
        # NB: TreeNode apparently needs unicode files
        with codecs.open( options.tree, 'rb', 'utf-8' ) as fh:
            extra_kwds['tree'] = TreeNode.read( fh )

    bd_dm = beta_diversity( options.distance_metric, counts, ids=sample_names, **extra_kwds )
    bd_dm.write( options.output )


if __name__ == "__main__":
    __main__()
