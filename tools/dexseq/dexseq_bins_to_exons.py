import argparse
import os
import re
import sys

'''
    File name:          dexseq_bins_to_exons.py
    Author:             Pavankumar Videm
    Date created:       25.07.2016
    Date last modified: 15.03.2017
    Purpose:            DEXSeq_prepare_annotation.py is used to prepare the DEXSeq compatible annotation
                        (flattened gtf file) from original GTF. In this process, the exons that appear multiple times,
                        once for each transcript are collapsed to so called "exon counting bins". Counting bins for
                        parts of exons arise when an exonic region appears with different boundaries in different
                        transcripts. The resulting flattened gtf file contains pseudo exon ids per gene instead of per
                        transcript. Assuming the  counting bins as actual exons leads to wrong conclusions. This script
                        annotates each DEXSeq result entry with geneid:transcriptid:exonnumber.
'''


def gtf_to_bed(gtf):
    """
    Parse the GTF file and write to a BED file with unique description formed by geneid, transcriptid and exon number.
    It assumed the exons are sorted by position for each gene.
    :param gtf: Initial GTF file used for DEXSeq prepare. Make sure that exonic features have an exon_number attribute.
    :return: none
    """
    fh_gtf = open(gtf, "r")
    fh_bed = open("annotation.bed", "w")
    prev_transcript = ""
    c = 0
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
        elif "transcript_id" not in d_description:
            sys.stderr.write("attribute \"transcript_id\" not found for gene :" + d_description["gene_id"])
            continue
        c += 1
        if prev_transcript != "" and prev_transcript != d_description["transcript_id"]:
            c = 1
        prev_transcript = d_description["transcript_id"]
        desc = d_description["gene_id"] + ":" + d_description["transcript_id"] + ":" + c
        fh_bed.write("\t".join([pos, desc, "0", strand]) + "\n")
    fh_gtf.close()
    fh_bed.close()


def main():
    parser = argparse.ArgumentParser(description='Maps the DEXSeq couting bins back to the transcript and exon ids')
    parser.add_argument('-i', '--input', help='DEXSeq output', required=True)
    parser.add_argument('-g', '--gtf', help='Original annotation gtf file used for DEXSeq prepare tool', required=True)
    parser.add_argument('-o', '--output', help='Output table', required=True)
    args = parser.parse_args()
    print("DEXSeq output file: %s" % args.input)
    print("Annotation file: %s" % args.gtf)
    print("Mapped output file: %s" % args.output)

    # create a BED entry for each line in the DEXSeq output file with the following entries and then sort it by
    # chromosome and start position
    # chr, start, end, geneids:dexseqbin, 0, strand
    os.system("awk '{print $12\"\t\"$13\"\t\"$14\"\t\"$1\"\t0\t\"$16}' "
              + args.input +
              " | sort -k1,1 -k2n,2 > output.bed")

    gtf_to_bed(args.gtf)

    # interset the DEXseq couting bins with exons in the GTF file
    # overlaped positions can be later used to infer which bin corresponds to which exon
    os.system("intersectBed -wo -s -a output.bed -b annotation.bed > overlap.txt")

    d_binexon = {}
    fh_overlap = open("overlap.txt", "r")
    for line in fh_overlap:
        binid = line.split('\t')[3]
        exonid = line.split('\t')[9]
        d_binexon.setdefault(binid, []).append(exonid)
    fh_overlap.close()

    fh_input = open(args.input, "r")
    fh_output = open(args.output, "w")
    for line in fh_input:
        binid = line.split('\t')[0]
        fh_output.write(line.rstrip('\n') + '\t' + ','.join(d_binexon[binid]) + '\n')
    fh_input.close()
    fh_output.close()

if __name__ == "__main__":
    main()
