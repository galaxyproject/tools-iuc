#!/usr/bin/env python3

"""
Converts a GFF produced by exonerate into a more standard GFF3 (e.g. usable in JBrowse)
"""

import argparse
import sys

from BCBio import GFF
from Bio.SeqFeature import FeatureLocation, SeqFeature

parser = argparse.ArgumentParser()
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('outfile', nargs='?', type=argparse.FileType('a'), default=sys.stdout)
args = parser.parse_args()


scaffs = []
gene_number = 0
for scaff in GFF.parse(args.infile):
    scaff.annotations = {}
    scaff.seq = ""
    kept_features = []
    current_gene = None
    exon_number = 0
    last_utr = None

    for feature in scaff.features:

        if feature.type == "gene":
            gene_number += 1
            mrna_feature = SeqFeature(FeatureLocation(feature.location.start, feature.location.end), type="mRNA", strand=feature.location.strand)
            mrna_feature.sub_features = []
            mrna_feature.qualifiers['source'] = feature.qualifiers['source']
            mrna_id = "mRNA_" + str(gene_number)
            mrna_feature.qualifiers['ID'] = mrna_id
            feature.sub_features = [mrna_feature]
            feature.qualifiers['ID'] = "gene_" + str(gene_number)
            if 'gene_orientation' in feature.qualifiers:
                del feature.qualifiers['gene_orientation']
            if current_gene:
                kept_features.append(current_gene)

            current_gene = feature
            exon_number = 0
            last_utr = None

        elif feature.type == 'utr5':
            feature.type = 'five_prime_UTR'
            feature.qualifiers['ID'] = '%s_five_prime_UTR' % (mrna_id)
            mrna_feature.sub_features.append(feature)
            last_utr = {'start': feature.location.start, 'end': feature.location.end}

        elif feature.type == 'utr3':
            feature.type = 'three_prime_UTR'
            feature.qualifiers['ID'] = '%s_three_prime_UTR' % (mrna_id)
            mrna_feature.sub_features.append(feature)
            last_utr = {'start': feature.location.start, 'end': feature.location.end}

        elif feature.type == 'exon':
            exon_number += 1
            feature.qualifiers['ID'] = '%s_exon_%s' % (mrna_id, exon_number)
            mrna_feature.sub_features.append(feature)

            if last_utr is None:
                cds_feature = SeqFeature(FeatureLocation(feature.location.start, feature.location.end), type="CDS", strand=feature.location.strand)
                cds_feature.sub_features = []
                cds_feature.qualifiers['source'] = feature.qualifiers['source']
                cds_feature.qualifiers['ID'] = mrna_id + "_CDS"
                mrna_feature.sub_features.append(cds_feature)
            elif feature.location.start != last_utr['start'] or feature.location.end != last_utr['end']:
                if feature.location.start > last_utr['start']:
                    cds_feature = SeqFeature(FeatureLocation(feature.location.start, last_utr['start']), type="CDS", strand=feature.location.strand)
                    cds_feature.sub_features = []
                    cds_feature.qualifiers['source'] = feature.qualifiers['source']
                    cds_feature.qualifiers['ID'] = mrna_id + "_CDS"
                    mrna_feature.sub_features.append(cds_feature)
                if feature.location.end < last_utr['end']:
                    cds_feature = SeqFeature(FeatureLocation(feature.location.end, last_utr['end']), type="CDS", strand=feature.location.strand)
                    cds_feature.sub_features = []
                    cds_feature.qualifiers['source'] = feature.qualifiers['source']
                    cds_feature.qualifiers['ID'] = mrna_id + "_CDS"
                    mrna_feature.sub_features.append(cds_feature)

            last_utr = None

        elif feature.type == 'similarity':
            if current_gene is None:
                # We haven't seen any gene, just convert similarity to match
                feature.type = 'match'
                kept_features.append(feature)

            last_utr = None

        elif feature.type not in ['splice3', 'splice5', 'similarity', 'intron']:
            mrna_feature.sub_features.append(feature)
            last_utr = None

    # For the last one
    if current_gene:
        kept_features.append(current_gene)

    scaff.features = kept_features

    if len(kept_features):
        GFF.write([scaff], args.outfile)
