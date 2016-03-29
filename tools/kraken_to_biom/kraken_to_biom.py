import argparse
import tempfile
import json
import os
import datetime
import subprocess


def main( args ):
    input_files = []
    for filename, label in zip( args.files.split( ',' ), args.labels.split( ',' ) ):
        fh, translated_filename = tempfile.mkstemp( suffix='.tab', prefix=filename )
        os.close( fh )
        with open( translated_filename, 'w' ) as fh:
            names_path = os.path.join( args.taxonomy_path, 'names.dmp' )
            nodes_path = os.path.join( args.taxonomy_path, 'nodes.dmp' )
            command = 'cat %s | taxonomy-reader %s %s 2' % ( filename, names_path, nodes_path )
            subprocess.call( command, stdout=fh, shell=True )
        input_files.append( ( filename, translated_filename, label ) )

    biom_data = dict( id=None,
                      format='Biological Observation Matrix 0.9.1-dev',
                      format_url='http://biom-format.org/documentation/format_versions/biom-1.0.html',
                      type='Taxon table',
                      generated_by='Galaxy',
                      date=datetime.datetime.now().isoformat(),
                      matrix_type='sparse',
                      matrix_element_type='int',
                      rows=[],
                      columns=[],
                      data=[] )

    tax_ids = {}
    taxa = []

    for original_filename, filename, label in input_files:
        biom_data['columns'].append( { 'id': label, 'metadata': None } )
        lines = file( filename, 'r' ).readlines()
        filtered_taxons = cluster_taxonomy( [ line.strip() for line in lines ], args.cluster )
        for line in lines:
            line_parts = line.split( '\t' )
            taxonomy = line_parts[ 5: ]
            tax_id = line_parts[ 2 ]
            if tax_id == '1':
                continue
            if filtered_taxons[ tax_id ] not in tax_ids:
                tax_ids[ filtered_taxons[ tax_id ] ] = { label: 1 }
            else:
                if label in tax_ids[ filtered_taxons[ tax_id ] ]:
                    tax_ids[ filtered_taxons[ tax_id ] ][ label ] += 1
                else:
                    tax_ids[ filtered_taxons[ tax_id ] ][ label ] = 1

    for idx, taxon in enumerate( tax_ids.keys() ):
        clean_taxonomy = []
        for tr in taxon.split(';'):
            if tr == 'n':
                continue
            clean_taxonomy.append( tr )
        if len( clean_taxonomy ) == 0:
            continue
        biom_data['rows'].append( { 'id': str( idx ), 'metadata': { 'taxonomy': clean_taxonomy, 'tax_id': taxon } } )

    rows = []

    for row in biom_data[ 'rows' ]:
        columns = []
        for column in biom_data[ 'columns' ]:
            tax_id = row[ 'metadata' ][ 'tax_id' ]
            label = column[ 'id' ]
            if tax_id in tax_ids:
                if label in tax_ids[ tax_id ]:
                    columns.append( tax_ids[ tax_id ][ label ] )
                else:
                    columns.append(0)
            else:
                columns.append(0)
        rows.append(columns)

    biom_data['shape'] = [ len( biom_data['rows'] ), len( biom_data['columns'] ) ]

    for i, row in enumerate( rows ):
        for j, column in enumerate( rows[ i ] ):
            if column != 0:
                biom_data[ 'data' ].append( [ i, j, column ] )

    file( args.output, 'w' ).write( json.dumps( biom_data, indent=2 ) )


def cluster_taxonomy( taxonomies, cluster_by ):
    taxonomy_mapping = {}
    taxonomy_filter = [ 'Classified', 'Read name', 'Tax id', 'Length', 'K-mers',
             'Superkingdom', 'Kingdom', 'Subkingdom', 'Superphylum', 'Phylum', 'Subphylum',
             'Superclass', 'Class', 'Subclass', 'Superorder', 'Order', 'Suborder',
             'Superfamily', 'Family', 'Subfamily', 'Tribe', 'Subtribe', 'Genus',
             'Subgenus', 'Species', 'Subspecies' ].index( cluster_by )
    for taxonomy in taxonomies:
        elements = taxonomy.split('\t')
        tax_id = elements[2]
        if len( elements ) == 0:
            continue
        if elements[0] == 'U':
            continue
        if tax_id == '1':
            continue
        tax_filter = ';'.join( elements[ 6:taxonomy_filter ] )
        if tax_id not in taxonomy_mapping:
            taxonomy_mapping[ tax_id ] = tax_filter
    return taxonomy_mapping


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    taxa = ( 'Superkingdom', 'Kingdom', 'Subkingdom', 'Superphylum', 'Phylum', 'Subphylum',
             'Superclass', 'Class', 'Subclass', 'Superorder', 'Order', 'Suborder',
             'Superfamily', 'Family', 'Subfamily', 'Tribe', 'Subtribe', 'Genus',
             'Subgenus', 'Species', 'Subspecies' )
    parser.add_argument( '--files', help="Comma-separated list of Kraken classifications", dest='files', action='store' )
    parser.add_argument( '--labels', help="Column labels", dest='labels', action='store' )
    parser.add_argument( '--taxonomy_path', help="Column labels", dest='taxonomy_path', action='store' )
    parser.add_argument( '--output', help="Output file", dest='output', action='store' )
    parser.add_argument( '--cluster', help="Combine entries below this taxonomic rank", dest='cluster', action='store', choices=taxa, default='Superclass' )
    args = parser.parse_args()
    exit(main(args))
