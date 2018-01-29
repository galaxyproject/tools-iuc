#!/usr/bin/env python

# Reports a summary of Kraken's results
# and optionally creates a newick Tree
# Copyright (c) 2016 Daniel Blankenberg
# Licensed under the Academic Free License version 3.0
# https://github.com/blankenberg/Kraken-Taxonomy-Report

from __future__ import print_function

import optparse
import os
import re
import sys

__VERSION__ = '0.0.2'

__URL__ = "https://github.com/blankenberg/Kraken-Taxonomy-Report"

# Rank names were pulled from ncbi nodes.dmp on 02/02/2016
# cat nodes.dmp | cut -f 5 | sort | uniq
# "root" is added manually
NO_RANK_NAME = "no rank"
RANK_NAMES = [ NO_RANK_NAME,
               "root",
               "superkingdom",
               "kingdom",
               "subkingdom",
               "superphylum",
               "phylum",
               "subphylum",
               "superclass",
               "class",
               "subclass",
               "infraclass",
               "superorder",
               "order",
               "suborder",
               "infraorder",
               "parvorder",
               "superfamily",
               "family",
               "subfamily",
               "tribe",
               "subtribe",
               "genus",
               "subgenus",
               "species group",
               "species subgroup",
               "species",
               "subspecies",
               "varietas",
               "forma" ]
# NB: We put 'no rank' at top of list for generating trees, due to e.g.
# root (root) -> cellular organisms (no rank) -> bacteria (superkingdom)

RANK_NAME_TO_INTS = dict( [ (y, x) for (x, y) in enumerate( RANK_NAMES ) ] )
RANK_NAMES_INTS = range( len( RANK_NAMES ) )

NO_RANK_INT = RANK_NAMES.index( NO_RANK_NAME )
NO_RANK_CODE = 'n'

PRIMARY_RANK_NAMES = [ 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom' ]
RANK_INT_TO_CODE = {}
for name in PRIMARY_RANK_NAMES:
    RANK_INT_TO_CODE[ RANK_NAMES.index( name ) ] = name[0]
RANK_INT_TO_CODE[ RANK_NAMES.index( 'superkingdom' ) ] = 'd'
PRIMARY_RANK_NAMES.append( 'superkingdom' )

NAME_STUB = "%s__%s"
NAME_RE = re.compile( "(\t| |\||\.;)" )
NAME_REPL = "_"


def get_kraken_db_path( db ):
    assert db, ValueError( "You must provide a kraken database" )
    k_db_path = os.getenv('KRAKEN_DB_PATH', None )
    if k_db_path:
        db = os.path.join( k_db_path, db )
    return db


def load_taxonomy( db_path, sanitize_names=False ):
    child_lists = {}
    name_map = {}
    rank_map = {}
    names = {}  # Store names here to look for duplicates (id, True/False name fixed)
    with open( os.path.join( db_path, "taxonomy/names.dmp" ) ) as fh:
        for line in fh:
            line = line.rstrip( "\n\r" )
            if line.endswith( "\t|" ):
                line = line[:-2]
            fields = line.split( "\t|\t" )
            node_id = fields[0]
            name = fields[1]
            if sanitize_names:
                name = NAME_RE.sub( NAME_REPL, name )
            name_type = fields[3]
            if name_type == "scientific name":
                if name in names:
                    print( 'Warning: name "%s" found at node "%s" but already exists originally for node "%s".' % ( name, node_id, names[name][0] ), file=sys.stderr )
                    new_name = "%s_%s" % ( name, node_id )
                    print( 'Transforming node "%s" named "%s" to "%s".' % ( node_id, name, new_name ), file=sys.stderr )
                    assert new_name not in names, 'Transformed Name "%s" already exists. Cannot recover at this time.' % new_name
                    if not names[name][1]:
                        orig_new_name = "%s_%s" % ( name, names[name][0] )
                        print( 'Transforming node "%s" named "%s" to "%s".' % ( names[name][0], name, orig_new_name ), file=sys.stderr )
                        assert orig_new_name not in names, 'Transformed Name "%s" already exists. Cannot recover at this time.' % orig_new_name
                        name_map[names[name][0]] = orig_new_name
                        names[name] = ( names[name][0], True )
                    name = new_name
                else:
                    names[name] = ( node_id, False )
                name_map[ node_id ] = name

    with open( os.path.join( db_path, "taxonomy/nodes.dmp" ) ) as fh:
        for line in fh:
            line = line.rstrip( "\n\r" )
            fields = line.split( "\t|\t" )
            node_id = fields[0]
            parent_id = fields[1]
            rank = RANK_NAME_TO_INTS.get( fields[2].lower(), None )
            if rank is None:
                # This should never happen, unless new taxonomy ranks are created
                print( 'Unrecognized rank: Node "%s" is "%s", setting to "%s"' % ( node_id, fields[2], NO_RANK_NAME ), file=sys.stderr )
                rank = NO_RANK_INT
            if node_id == '1':
                parent_id = '0'
            if parent_id not in child_lists:
                child_lists[ parent_id ] = []
            child_lists[ parent_id ].append( node_id )
            rank_map[node_id] = rank
    return ( child_lists, name_map, rank_map )


def dfs_summation( node, counts, child_lists ):
    children = child_lists.get( node, None )
    if children:
        for child in children:
            dfs_summation( child, counts, child_lists )
            counts[ node ] = counts.get( node, 0 ) + counts.get( child, 0 )


def dfs_report( node, file_data, hit_taxa, rank_map, name_map, child_lists, output_lines, options, name=None, tax=None ):
    rank_int = rank_map[node]
    code = RANK_INT_TO_CODE.get( rank_int, NO_RANK_CODE )
    if ( code != NO_RANK_CODE or options.intermediate ) and ( options.show_zeros or node in hit_taxa):
        if name is None:
            name = ""
        else:
            name = "%s|" % name
        if tax is None:
            tax = ''
        else:
            tax = "%s;" % tax
        sanitized_name = name_map[ node ]
        name_stub = NAME_STUB % ( code, sanitized_name )
        name = name + name_stub
        tax = tax + name_stub
        if options.name_id:
            output = node
        elif options.name_long:
            output = name
        else:
            output = sanitized_name
        for val in file_data:
            output = "%s\t%i" % ( output, val.get( node, 0 ) )
        if options.show_rank:
            output = "%s\t%s" % ( output, RANK_NAMES[ rank_int ] )
        if options.taxonomy:
            output = "%s\t%s" % ( output, tax )
        output_lines[ rank_int ].append( output )
    children = child_lists.get( node )
    if children:
        for child in children:
            dfs_report( child, file_data, hit_taxa, rank_map, name_map, child_lists, output_lines, options, name=name, tax=tax )


def write_tree( child_lists, name_map, rank_map, options, branch_length=1 ):
    # Uses Biopython, only load if making tree
    import Bio.Phylo
    from Bio.Phylo import BaseTree

    def _get_name( node_id ):
        if options.name_id:
            return node_id
        return name_map[node_id]
    nodes = {}
    root_node_id = child_lists["0"][0]
    nodes[root_node_id] = BaseTree.Clade( name=_get_name( root_node_id), branch_length=branch_length )

    def recurse_children( parent_id ):
        if options.cluster is not None and rank_map[parent_id] == options.cluster:
            # Short circuit if we found our rank, prevents 'hanging' no ranks from being output
            # e.g. clustering by "species" (Escherichia coli), but have "no rank" below (Escherichia coli K-12) in test_db
            return
        if parent_id not in nodes:
            nodes[parent_id] = BaseTree.Clade( name=_get_name( parent_id ), branch_length=branch_length )
        for child_id in child_lists.get( parent_id, [] ):
            if options.cluster is None or ( rank_map[child_id] <= options.cluster  ):
                if child_id not in nodes:
                    nodes[child_id] = BaseTree.Clade(name=_get_name( child_id ), branch_length=branch_length)
                nodes[parent_id].clades.append(nodes[child_id])
                recurse_children( child_id )
    recurse_children( root_node_id )
    tree = BaseTree.Tree(root=nodes[root_node_id])
    Bio.Phylo.write( [tree], options.output_tree, 'newick' )


def __main__():
    parser = optparse.OptionParser( usage="%prog [options] file1 file...fileN" )
    parser.add_option( '-v', '--version', dest='version', action='store_true', default=False, help='print version and exit' )
    parser.add_option( '', '--show-zeros', dest='show_zeros', action='store_true', default=False, help='Show empty nodes' )
    parser.add_option( '', '--header-line', dest='header_line', action='store_true', default=False, help='Provide a header on output' )
    parser.add_option( '', '--intermediate', dest='intermediate', action='store_true', default=False, help='Intermediate Ranks' )
    parser.add_option( '', '--name-id', dest='name_id', action='store_true', default=False, help='Use Taxa ID instead of Name' )
    parser.add_option( '', '--name-long', dest='name_long', action='store_true', default=False, help='Use Long taxa ID instead of base name' )
    parser.add_option( '', '--taxonomy', dest='taxonomy', action='store_true', default=False, help='Output taxonomy in last column' )
    parser.add_option( '', '--cluster', dest='cluster', action='store', type="string", default=None, help='Cluster counts to specified rank' )
    parser.add_option( '', '--summation', dest='summation', action='store_true', default=False, help='Add summation of child counts to each taxa' )
    parser.add_option( '', '--sanitize-names', dest='sanitize_names', action='store_true', default=False, help='Replace special chars (\t| |\||\.;) with underscore (_)' )
    parser.add_option( '', '--show-rank', dest='show_rank', action='store_true', default=False, help='Output column with Rank name' )
    parser.add_option( '', '--db', dest='db', action='store', type="string", default=None, help='Name of Kraken database' )
    parser.add_option( '', '--output', dest='output', action='store', type="string", default=None, help='Name of output file' )
    parser.add_option( '', '--output-tree', dest='output_tree', action='store', type="string", default=None, help='Name of output file to place newick tree' )
    (options, args) = parser.parse_args()
    if options.version:
        print( "Kraken Taxonomy Report (%s) version %s" % ( __URL__, __VERSION__ ), file=sys.stderr )
        sys.exit()
    if not args:
        print( parser.get_usage(), file=sys.stderr )
        sys.exit()

    if options.cluster:
        cluster_name = options.cluster.lower()
        cluster = RANK_NAME_TO_INTS.get( cluster_name, None )
        assert cluster is not None, ValueError( '"%s" is not a valid rank for clustering.' % options.cluster )
        if cluster_name not in PRIMARY_RANK_NAMES:
            assert options.intermediate, ValueError( 'You cannot cluster by "%s", unless you enable intermediate ranks.' % options.cluster )
        ranks_to_report = [ cluster ]
        options.cluster = cluster
        # When clustering we need to do summatation
        options.summation = True
    else:
        options.cluster = None  # make empty string into None
        ranks_to_report = RANK_NAMES_INTS

    if options.output:
        output_fh = open(options.output, 'w')
    else:
        output_fh = sys.stdout

    db_path = get_kraken_db_path( options.db )
    ( child_lists, name_map, rank_map ) = load_taxonomy( db_path, sanitize_names=options.sanitize_names )
    file_data = []
    hit_taxa = []
    for input_filename in args:
        taxo_counts = {}
        with open( input_filename ) as fh:
            for line in fh:
                fields = line.split( "\t" )
                taxo_counts[ fields[2] ] = taxo_counts.get( fields[2], 0 ) + 1
        clade_counts = taxo_counts.copy()  # fixme remove copying?
        if options.summation:
            dfs_summation( '1', clade_counts, child_lists )
        for key, value in clade_counts.items():
            if value and key not in hit_taxa:
                hit_taxa.append( key )
        file_data.append( clade_counts )

    if options.header_line:
        output_fh.write( "#ID\t" )
        output_fh.write( "\t".join( args ) )
        if options.show_rank:
            output_fh.write( "\trank" )
        if options.taxonomy:
            output_fh.write( "\ttaxonomy" )
        output_fh.write( '\n' )

    output_lines = dict( [ ( x, [] ) for x in RANK_NAMES_INTS ] )
    dfs_report( '1', file_data, hit_taxa, rank_map, name_map, child_lists, output_lines, options, name=None, tax=None )

    for rank_int in ranks_to_report:
        for line in output_lines.get( rank_int, [] ):
            output_fh.write( line )
            output_fh.write( '\n' )
    fh.close()
    if options.output_tree:
        write_tree( child_lists, name_map, rank_map, options )


if __name__ == "__main__":
    __main__()
