import os
import sys
import argparse
from collections import defaultdict
from BCBio import GFF


def strandardize(strand):
    if str(strand) == '-1':
        strand = '-'
    elif str(strand) == '1':
        strand = '+'
    return strand


def gtf_to_dict(f_gff, feat_type, idattr, attributes):
    """
    It reads only exonic features because not all GFF files contain gene and trascript features. From the exonic
    features it extracts gene names, biotypes, start and end positions. If any of these attributes do not exit
    then they are set to NA.
    """
    annotation = defaultdict(lambda: defaultdict(lambda: 'NA'))
    exon_pos = defaultdict(lambda: defaultdict(dict))

    with open(f_gff) as gff_handle:
        for rec in GFF.parse(gff_handle, limit_info=dict(gff_type=[feat_type]), target_lines=1):
            for sub_feature in rec.features:
                start = sub_feature.location.start
                end = sub_feature.location.end
                strand = strandardize(sub_feature.location.strand)
                try:
                    geneid = sub_feature.qualifiers[idattr][0]
                except KeyError:
                    print("No '" + idattr + "' attribute found for the feature at position "
                          + rec.id + ":" + str(start) + ":" + str(end) + ". Please check your GTF/GFF file.")
                    continue

                annotation[geneid]['chr'] = rec.id
                annotation[geneid]['strand'] = strand
                if annotation[geneid]['start'] == 'NA' or start <= int(annotation[geneid]['start']):
                    annotation[geneid]['start'] = start
                if annotation[geneid]['end'] == 'NA' or end >= int(annotation[geneid]['end']):
                    annotation[geneid]['end'] = end

                for attr in attributes:
                    if attr in annotation[geneid]:
                        continue
                    try:
                        annotation[geneid][attr] = sub_feature.qualifiers[attr][0]
                    except KeyError:
                        annotation[geneid][attr] = 'NA'
    gff_handle.close()

    return annotation


def main():
    parser = argparse.ArgumentParser(description='Annotate DESeq2/DEXSeq tables with more information from GTF files')
    parser.add_argument('-in', '--input', required=True, help='DESeq2 output', )
    parser.add_argument('-g', '--gtf', required=True, help='Annotation GTF file used for couting')
    parser.add_argument('-t', '--type', default='exon', required=False,
                        help='feature type (3rd column in GFF file) to be used (default: exon)')
    parser.add_argument('-i', '--idattr', default='gene_id', required=False,
                        help='GFF attribute to be used as feature ID. '
                             'This should match the first column of DESeq2 output(default: geneid)')
    parser.add_argument('-a', '--attributes', default = 'gene_biotype, gene_name', required=False,
                        help='Comma separated attributes to include in output')
    parser.add_argument('-o', '--output', required=True, help='Output file path')
    args = parser.parse_args()

    print("DE(X)Seq output file     : %s" % args.input)
    print("Annotation file          : %s" % args.gtf)
    print("Feature type             : %s" % args.type)
    print("ID attribute             : %s" % args.idattr)
    print("Attributes to include    : %s" % args.attributes)
    print("Annotated output file    : %s" % args.output)

    attr = [x.strip() for x in args.attributes.split(',')]

    annotation = gtf_to_dict(args.gtf, args.type, args.idattr, attr)

    fh_input = open(args.input, "r")
    fh_output = open(args.output, "w")
    for line in fh_input:
        # Append the extra information from GTF to DESeq2 output
        geneid = line.split('\t')[0]
        annot = [str(annotation[geneid]['chr']),
                 str(annotation[geneid]['start']),
                 str(annotation[geneid]['end']),
                 str(annotation[geneid]['strand'])]
        for a in attr:
            annot.append(annotation[geneid][a])
        fh_output.write(line.rstrip('\n') + '\t' + '\t'.join(annot) + '\n')
    fh_input.close()
    fh_output.close()


if __name__ == "__main__":
    main()
