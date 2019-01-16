import argparse
import os
from collections import defaultdict

from BCBio import GFF


def strandardize(strand):
    if str(strand) == '-1':
        strand = '-'
    elif str(strand) == '1':
        strand = '+'
    return strand


def gff_to_dict(f_gff, feat_type, idattr, txattr, attributes, input_type):
    """
    It reads only exonic features because not all GFF files contain gene and trascript features. From the exonic
    features it extracts gene names, biotypes, start and end positions. If any of these attributes do not exit
    then they are set to NA.
    """
    annotation = defaultdict(lambda: defaultdict(lambda: 'NA'))
    exon_pos = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    tx_info = defaultdict(lambda: defaultdict(str))

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
                # extract exon information only in case of dexseq output
                if input_type != "dexseq":
                    continue
                try:
                    txid = sub_feature.qualifiers[txattr][0]
                    tx_info[txid]['chr'] = rec.id
                    tx_info[txid]['strand'] = strand
                    exon_pos[txid][int(start)][int(end)] = 1
                except KeyError:
                    print("No '" + txattr + "' attribute found for the feature at position " + rec.id + ":" + str(
                        start) + ":" + str(end) + ". Please check your GTF/GFF file.")
                    pass

    bed_entries = []
    # create BED lines only for deseq output
    if input_type == "dexseq":
        for txid in exon_pos.keys():
            starts = sorted(exon_pos[txid])
            strand = tx_info[txid]['strand']
            if strand == '-':
                starts = reversed(starts)
            for c, start in enumerate(starts, 1):
                ends = sorted(exon_pos[txid][start])
                if strand == '-':
                    ends = reversed(ends)
                for end in ends:
                    bed_entries.append('\t'.join([tx_info[txid]['chr'], str(start), str(end),
                                                  txid + ':' + str(c), '0', strand]))

    return annotation, bed_entries


def main():
    parser = argparse.ArgumentParser(description='Annotate DESeq2/DEXSeq tables with information from GFF/GTF files')
    parser.add_argument('-in', '--input', required=True,
                        help='DESeq2/DEXSeq output. It is allowed to have extra information, '
                             'but make sure that the original output columns are not altered')
    parser.add_argument('-m', '--mode', required=True, choices=["deseq2", "dexseq"], default='deseq2',
                        help='Input file type')
    parser.add_argument('-g', '--gff', required=True, help='The same annotation GFF/GTF file used for couting')
    parser.add_argument('-t', '--type', default='exon', required=False,
                        help='feature type (3rd column in GFF file) to be used (default: exon)')
    parser.add_argument('-i', '--idattr', default='gene_id', required=False,
                        help='GFF attribute to be used as feature ID. '
                             'This should match the first column of DESeq2 output(default: geneid)')
    parser.add_argument('-x', '--txattr', default='transcript_id', required=False,
                        help='GFF attribute to be used as transcript ID. Used for DEXSeq output only.'
                             'This should match the first column of DESeq2 output(default: transcript_id)')
    parser.add_argument('-a', '--attributes', default='gene_biotype, gene_name', required=False,
                        help='Comma separated attributes to include in output. Default: gene_biotype, gene_name')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    args = parser.parse_args()

    print("DE(X)Seq output file     : %s" % args.input)
    print("Input file type          : %s" % args.mode)
    print("Annotation file          : %s" % args.gff)
    print("Feature type             : %s" % args.type)
    print("ID attribute             : %s" % args.idattr)
    print("Transcript attribute     : %s" % args.txattr)
    print("Attributes to include    : %s" % args.attributes)
    print("Annotated output file    : %s" % args.output)

    attr = [x.strip() for x in args.attributes.split(',')]
    annotation, bed_entries = gff_to_dict(args.gff, args.type, args.idattr, args.txattr, attr, args.mode)

    d_binexon = {}
    skip_exon_annotation = False

    if args.mode == "dexseq":
        with open(args.input) as fh_input, open("input.bed", "w") as fh_input_bed:
            for line in fh_input:
                f = line.split('\t')
                fh_input_bed.write('\t'.join([f[11], f[12], f[13], f[0], "0", f[15]]) + "\n")

        if len(bed_entries) == 0 and args.mode == "dexseq":
            print("It seems there are no transcript ids present in GFF file. Skipping exon annotation.")
            skip_exon_annotation = True

        if not skip_exon_annotation:
            with open("annotation.bed", "w") as fh_annotation_bed:
                for line in bed_entries:
                    fh_annotation_bed.write(line + "\n")

            # interset the DEXseq couting bins with exons in the GFF file
            # overlaped positions can be later used to infer which bin corresponds to which exon
            os.system("intersectBed -wo -s -a input.bed -b annotation.bed > overlap.txt")

            with open("overlap.txt") as fh_overlap:
                for line in fh_overlap:
                    binid = line.split('\t')[3]
                    exonid = line.split('\t')[9]
                    d_binexon.setdefault(binid, []).append(exonid)

    with open(args.input) as fh_input, open(args.output, 'w') as fh_output:
        for line in fh_input:
            annot = []
            # Append the extra information from GFF to DESeq2 output
            if args.mode == "deseq2":
                geneid = line.split('\t')[0]
                annot = [str(annotation[geneid]['chr']),
                         str(annotation[geneid]['start']),
                         str(annotation[geneid]['end']),
                         str(annotation[geneid]['strand'])]
                for a in attr:
                    annot.append(annotation[geneid][a])
            # DEXSeq exonic bins might originate from aggrigating multiple genes. They are are separated by '+'
            # Append the attributes from the GFF but keep the order of the aggregated genes and use '+'
            # Aappend the transcript id and exon number from the annotation that correspond to the DEXseq counting bins
            elif args.mode == "dexseq":
                geneids = line.split('\t')[1].split('+')
                for a in attr:
                    tmp = []
                    for geneid in geneids:
                        tmp.append(str(annotation[geneid][a]))
                    annot.append('+'.join(tmp))
                if not skip_exon_annotation:
                    binid = line.split('\t')[0]
                    try:
                        annot.append(','.join(sorted(set(d_binexon[binid]))))
                    except KeyError:
                        annot.append('NA')
            fh_output.write(line.rstrip('\n') + '\t' + '\t'.join(annot) + '\n')


if __name__ == "__main__":
    main()
