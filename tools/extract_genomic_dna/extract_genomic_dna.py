#!/usr/bin/env python
from __future__ import print_function

import argparse
import os

import bx.seq.nib
import bx.seq.twobit
from bx.intervals.io import Comment, Header

import extract_genomic_dna_utils as egdu  # noqa: I100,I202

parser = argparse.ArgumentParser()
parser.add_argument('--input_format', dest='input_format', help="Input dataset format")
parser.add_argument('--input', dest='input', help="Input dataset")
parser.add_argument('--genome', dest='genome', help="Input dataset genome build")
parser.add_argument('--interpret_features', dest='interpret_features', default=None, help="Interpret features if input format is gff")
parser.add_argument('--columns', dest='columns', help="Columns to use in input file")
parser.add_argument('--reference_genome_source', dest='reference_genome_source', help="Source of reference genome file")
parser.add_argument('--reference_genome', dest='reference_genome', help="Reference genome file")
parser.add_argument('--output_format', dest='output_format', help="Output format")
parser.add_argument('--fasta_header_type', dest='fasta_header_type', default=None, help="Fasta header format")
parser.add_argument('--fasta_header_delimiter', dest='fasta_header_delimiter', default=None, help="Fasta header field delimiter")
parser.add_argument('--output', dest='output', help="Output dataset")
args = parser.parse_args()

input_is_gff = args.input_format == 'gff'
interpret_features = input_is_gff and args.interpret_features == "yes"
if len(args.columns.split(',')) == 5:
    # Bed file.
    chrom_col, start_col, end_col, strand_col, name_col = egdu.parse_cols_arg(args.columns)
else:
    # Gff file.
    chrom_col, start_col, end_col, strand_col = egdu.parse_cols_arg(args.columns)
    name_col = False

if args.reference_genome_source == "history":
    seq_path = egdu.convert_to_twobit(args.reference_genome)
else:
    seq_path = args.reference_genome
seq_dir = os.path.split(seq_path)[0]

includes_strand_col = strand_col >= 0
strand = None
nibs = {}
skipped_lines = 0
first_invalid_line = 0
invalid_lines = []
warnings = []
warning = ''
twobitfile = None
line_count = 1
file_iterator = open(args.input)
if interpret_features:
    file_iterator = egdu.GFFReaderWrapper(file_iterator, fix_strand=False)
out = open(args.output, 'wt')

for feature in file_iterator:
    # Ignore comments, headers.
    if isinstance(feature, (Header, Comment)):
        line_count += 1
        continue
    name = ""
    if interpret_features:
        # Processing features.
        egdu.convert_gff_coords_to_bed(feature)
        chrom = feature.chrom
        start = feature.start
        end = feature.end
        strand = feature.strand
    else:
        # Processing lines, either interval or GFF format.
        line = feature.rstrip('\r\n')
        if line and not line.startswith("#"):
            fields = line.split('\t')
            try:
                chrom = fields[chrom_col]
                start = int(fields[start_col])
                end = int(fields[end_col])
                if name_col:
                    name = fields[name_col]
                if input_is_gff:
                    start, end = egdu.convert_gff_coords_to_bed([start, end])
                if includes_strand_col:
                    strand = fields[strand_col]
            except Exception:
                warning = "Invalid chrom, start or end column values. "
                warnings.append(warning)
                if not invalid_lines:
                    invalid_lines = egdu.get_lines(feature)
                    first_invalid_line = line_count
                skipped_lines += len(invalid_lines)
                continue
            if start > end:
                warning = "Invalid interval, start '%d' > end '%d'.  " % (start, end)
                warnings.append(warning)
                if not invalid_lines:
                    invalid_lines = egdu.get_lines(feature)
                    first_invalid_line = line_count
                skipped_lines += len(invalid_lines)
                continue
            if strand not in ['+', '-']:
                strand = '+'
            sequence = ''
        else:
            continue
    # Open sequence file and get sequence for feature/interval.
    if os.path.exists("%s/%s.nib" % (seq_dir, chrom)):
        if chrom in nibs:
            nib = nibs[chrom]
        else:
            nibs[chrom] = nib = bx.seq.nib.NibFile(open("%s/%s.nib" % (seq_path, chrom)))
        try:
            sequence = nib.get(start, end - start)
        except Exception:
            warning = "Unable to fetch the sequence from '%d' to '%d' for build '%s'. " % (start, end - start, args.genome)
            warnings.append(warning)
            if not invalid_lines:
                invalid_lines = egdu.get_lines(feature)
                first_invalid_line = line_count
            skipped_lines += len(invalid_lines)
            continue
    elif os.path.isfile(seq_path):
        if not(twobitfile):
            twobitfile = bx.seq.twobit.TwoBitFile(open(seq_path))
        try:
            if interpret_features:
                # Create sequence from intervals within a feature.
                sequence = ''
                for interval in feature.intervals:
                    sequence += twobitfile[interval.chrom][interval.start:interval.end]
            else:
                sequence = twobitfile[chrom][start:end]
        except Exception:
            warning = "Unable to fetch the sequence from '%d' to '%d' for chrom '%s'. " % (start, end - start, chrom)
            warnings.append(warning)
            if not invalid_lines:
                invalid_lines = egdu.get_lines(feature)
                first_invalid_line = line_count
            skipped_lines += len(invalid_lines)
            continue
    else:
        warning = "Chromosome by name '%s' was not found for build '%s'. " % (chrom, args.genome)
        warnings.append(warning)
        if not invalid_lines:
            invalid_lines = egdu.get_lines(feature)
            first_invalid_line = line_count
        skipped_lines += len(invalid_lines)
        continue
    if sequence == '':
        warning = "Chrom: '%s', start: '%d', end: '%d' is either invalid or not present in build '%s'. " % (chrom, start, end, args.genome)
        warnings.append(warning)
        if not invalid_lines:
            invalid_lines = egdu.get_lines(feature)
            first_invalid_line = line_count
        skipped_lines += len(invalid_lines)
        continue
    if includes_strand_col and strand == "-":
        sequence = egdu.reverse_complement(sequence)
    if args.output_format == "fasta":
        if input_is_gff:
            start, end = egdu.convert_bed_coords_to_gff([start, end])
        if args.fasta_header_type == "bedtools_getfasta_default":
            out.write(">%s\n" % egdu.get_bedtools_getfasta_default_header(str(chrom),
                                                                          str(start),
                                                                          str(end),
                                                                          strand,
                                                                          includes_strand_col))
        else:
            # args.fasta_header_type == "char_delimited":
            fields = [args.genome, str(chrom), str(start), str(end), strand]
            field_delimiter = egdu.get_fasta_header_delimiter(args.fasta_header_delimiter)
            meta_data = field_delimiter.join(fields)
            if name.strip():
                out.write(">%s %s\n" % (meta_data, name))
            else:
                out.write(">%s\n" % meta_data)
        c = 0
        sequence_length = len(sequence)
        while c < sequence_length:
            b = min(c + 50, sequence_length)
            out.write("%s\n" % str(sequence[c:b]))
            c = b
    else:
        # output_format == "interval".
        if interpret_features:
            meta_data = "\t".join([feature.chrom,
                                   "galaxy_extract_genomic_dna",
                                   "interval",
                                   str(feature.start),
                                   str(feature.end),
                                   feature.score,
                                   feature.strand,
                                   ".",
                                   egdu.gff_attributes_to_str(feature.attributes, "GTF")])
        else:
            # Here fields was set up around line 73.
            meta_data = "\t".join(fields)
        if input_is_gff:
            format_str = "%s seq \"%s\";\n"
        else:
            format_str = "%s\t%s\n"
        out.write(format_str % (meta_data, str(sequence)))
    # Update line count.
    if isinstance(feature, egdu.GFFFeature):
        line_count += len(feature.intervals)
    else:
        line_count += 1
out.close()

if warnings:
    warn_msg = "%d warnings, 1st is: " % len(warnings)
    warn_msg += warnings[0]
    print(warn_msg)
if skipped_lines:
    # Error message includes up to the first 10 skipped lines.
    print('Skipped %d invalid lines, 1st is #%d, "%s"' % (skipped_lines, first_invalid_line, '\n'.join(invalid_lines[:10])))

if args.reference_genome_source == "history":
    os.remove(seq_path)
