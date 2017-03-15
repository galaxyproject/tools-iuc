import argparse
import os
import re

'''
    File name:          deseq2_annotate.py
    Author:             Pavankumar Videm
    Date created:       31.07.2016
    Date last modified: 15.03.2017
    Purpose:            To add some extra useful information like gene names, biotype etc. from the GTF file to
                        DESeq2/DEXSeq output file.
'''


def gtf_to_dict(gtf, input_type):
    """
    It reads only exonic features because not all GTF files contain gene and trascript features. From the exonic
    features it extracts gene names, biotypes, start and end positions. If any of these attributes do not exit
    then they are set to NA.
    :param gtf: GTF file with exonic features
    :param input_type: Type of the input, string deseq2 or dexseq. Depending on the type different information is
    stored in the dictionary
    :return: dictionary containing geneid -> information. DEXSeq output already contains position information.
    Hence position is ignored for DEXSeq dictionary.
    """
    fh_gtf = open(gtf, "r")
    d_gene_startpos = {}
    d_gene_endpos = {}
    d_annotation = {}
    for line in fh_gtf:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if f[2] != "exon":
            continue
        pos = '\t'.join([f[0], f[3], f[4]])
        strand = f[6]
        d_description = dict(item.replace("\";", "").split(" \"") for item in filter(None, f[8].rstrip().split("\"; ")))
        if "gene_id" not in d_description:
            sys.stderr.write("attribute \"gene_id\" not found for position :" + pos)
            continue
        if d_description["gene_id"] in d_gene_startpos:
            if int(f[3]) < int(d_gene_startpos[d_description["gene_id"]]):
                d_gene_startpos[d_description["gene_id"]] = int(f[3])
            if int(f[4]) > int(d_gene_endpos[d_description["gene_id"]]):
                d_gene_endpos[d_description["gene_id"]] = int(f[4])
        else:
            d_gene_startpos[d_description["gene_id"]] = int(f[3])
            d_gene_endpos[d_description["gene_id"]] = int(f[4])
        if "gene_biotype" not in d_description:
            d_description["gene_biotype"] = "NA"
        if "gene_name" not in d_description:
            d_description["gene_name"] = "NA"
        desc = ""
        if input_type == "deseq2":
            desc = "\t".join([d_description["gene_biotype"],
                              d_description["gene_name"],
                              f[0],
                              str(d_gene_startpos[d_description["gene_id"]]),
                              str(d_gene_endpos[d_description["gene_id"]]),
                              strand])
        elif input_type == "dexseq":
            desc = d_description["gene_biotype"] + "\t" + d_description["gene_name"]
        d_annotation[d_description["gene_id"]] = desc
    fh_gtf.close()
    return d_annotation


def main():
    parser = argparse.ArgumentParser(description='Annotate DESeq2/DEXSeq tables with more information from GTF files')
    parser.add_argument('-i', '--input', help='DESeq2/DEXSeq output', required=True)
    parser.add_argument('-g', '--gtf', help='Original annotation GTF file used', required=True)
    parser.add_argument('-o', '--output', help='Output table', required=True)
    args = parser.parse_args()
    print("DE(X)Seq output file: %s" % args.input)
    print("Annotation file: %s" % args.gtf)
    print("Annotated output file: %s" % args.output)

    # DEseq2 output file contains 7 fields and DEXseq output file 27 fields. This inforamtion may change
    # in future versions of the tools.
    with open(args.input, "r") as f:
        line = f.readline()
        if len(line.split('\t')) == 7:
            input_type = "deseq2"
        elif len(line.split('\t')) == 27:
            input_type = "dexseq"
        else:
            sys.exit("Unkown file format! Please provide the unchanged DESeq2 or DEXSeq output.")

    d_annotation = gtf_to_dict(args.gtf, input_type)

    fh_input = open(args.input, "r")
    fh_output = open(args.output, "w")
    for line in fh_input:
        # Append the extra information from GTF to DESeq2 output
        if input_type == "deseq2":
            geneid = line.split('\t')[0]
            fh_output.write(line.rstrip('\n') + '\t' + d_annotation[geneid] + '\n')
        # DEXSeq exonic bins might originate from aggrigating multiple genes. They are are separated by '+'
        # Append the gene name and biotype from the GTF but keep the order of the aggregated genes and use '+'
        # concatenate the aggregated genes information
        elif input_type == "dexseq":
            l_geneids = line.split('\t')[1].split('+')
            l_genenames = []
            l_genetypes = []
            for geneid in l_geneids:
                [genename, genetype] = d_annotation[geneid].split('\t')
                l_genenames.append(genename)
                l_genetypes.append(genetype)
            fh_output.write(line.rstrip('\n') + '\t' + '+'.join(l_genetypes) + '\t' + '+'.join(l_genenames) + '\n')
    fh_input.close()
    fh_output.close()

if __name__ == "__main__":
    main()
