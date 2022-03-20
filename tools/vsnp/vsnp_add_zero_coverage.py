#!/usr/bin/env python

import argparse
import os
import re
import shutil

import pandas
import pysam
from Bio import SeqIO


def get_sample_name(file_path):
    base_file_name = os.path.basename(file_path)
    if base_file_name.find(".") > 0:
        # Eliminate the extension.
        return os.path.splitext(base_file_name)[0]
    return base_file_name


def get_coverage_df(bam_file):
    # Create a coverage dictionary.
    coverage_dict = {}
    coverage_list = pysam.depth(bam_file, split_lines=True)
    for line in coverage_list:
        chrom, position, depth = line.split('\t')
        coverage_dict["%s-%s" % (chrom, position)] = depth
    # Convert it to a data frame.
    coverage_df = pandas.DataFrame.from_dict(coverage_dict, orient='index', columns=["depth"])
    return coverage_df


def get_zero_df(reference):
    # Create a zero coverage dictionary.
    zero_dict = {}
    for record in SeqIO.parse(reference, "fasta"):
        chrom = record.id
        total_len = len(record.seq)
        for pos in list(range(1, total_len + 1)):
            zero_dict["%s-%s" % (str(chrom), str(pos))] = 0
    # Convert it to a data frame with depth_x
    # and depth_y columns - index is NaN.
    zero_df = pandas.DataFrame.from_dict(zero_dict, orient='index', columns=["depth"])
    return zero_df


def output_zc_vcf_file(base_file_name, vcf_file, zero_df, total_zero_coverage, output_vcf):
    column_names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample"]
    vcf_df = pandas.read_csv(vcf_file, sep='\t', header=None, names=column_names, comment='#')
    good_snp_count = len(vcf_df[(vcf_df['ALT'].str.len() == 1) & (vcf_df['REF'].str.len() == 1) & (vcf_df['QUAL'] > 150)])
    if total_zero_coverage > 0:
        header_file = "%s_header.csv" % base_file_name
        with open(header_file, 'w') as outfile:
            with open(vcf_file) as infile:
                for line in infile:
                    if re.search('^#', line):
                        outfile.write("%s" % line)
        vcf_df_snp = vcf_df[vcf_df['REF'].str.len() == 1]
        vcf_df_snp = vcf_df_snp[vcf_df_snp['ALT'].str.len() == 1]
        vcf_df_snp['ABS_VALUE'] = vcf_df_snp['CHROM'].map(str) + "-" + vcf_df_snp['POS'].map(str)
        vcf_df_snp = vcf_df_snp.set_index('ABS_VALUE')
        cat_df = pandas.concat([vcf_df_snp, zero_df], axis=1, sort=False)
        cat_df = cat_df.drop(columns=['CHROM', 'POS', 'depth'])
        cat_df[['ID', 'ALT', 'QUAL', 'FILTER', 'INFO']] = cat_df[['ID', 'ALT', 'QUAL', 'FILTER', 'INFO']].fillna('.')
        cat_df['REF'] = cat_df['REF'].fillna('N')
        cat_df['FORMAT'] = cat_df['FORMAT'].fillna('GT')
        cat_df['Sample'] = cat_df['Sample'].fillna('./.')
        cat_df['temp'] = cat_df.index.str.rsplit('-', n=1)
        cat_df[['CHROM', 'POS']] = pandas.DataFrame(cat_df.temp.values.tolist(), index=cat_df.index)
        cat_df = cat_df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'Sample']]
        cat_df['POS'] = cat_df['POS'].astype(int)
        cat_df = cat_df.sort_values(['CHROM', 'POS'])
        body_file = "%s_body.csv" % base_file_name
        cat_df.to_csv(body_file, sep='\t', header=False, index=False)
        with open(output_vcf, "w") as outfile:
            for cf in [header_file, body_file]:
                with open(cf, "r") as infile:
                    for line in infile:
                        outfile.write("%s" % line)
    else:
        shutil.move(vcf_file, output_vcf)
    return good_snp_count


def output_metrics_file(base_file_name, average_coverage, genome_coverage, good_snp_count, output_metrics):
    bam_metrics = [base_file_name, "", "%4f" % average_coverage, genome_coverage]
    vcf_metrics = [base_file_name, str(good_snp_count), "", ""]
    metrics_columns = ["File", "Number of Good SNPs", "Average Coverage", "Genome Coverage"]
    with open(output_metrics, "w") as fh:
        fh.write("# %s\n" % "\t".join(metrics_columns))
        fh.write("%s\n" % "\t".join(bam_metrics))
        fh.write("%s\n" % "\t".join(vcf_metrics))


def output_files(vcf_file, total_zero_coverage, zero_df, output_vcf, average_coverage, genome_coverage, output_metrics):
    base_file_name = get_sample_name(vcf_file)
    good_snp_count = output_zc_vcf_file(base_file_name, vcf_file, zero_df, total_zero_coverage, output_vcf)
    output_metrics_file(base_file_name, average_coverage, genome_coverage, good_snp_count, output_metrics)


def get_coverage_and_snp_count(bam_file, vcf_file, reference, output_metrics, output_vcf):
    coverage_df = get_coverage_df(bam_file)
    zero_df = get_zero_df(reference)
    coverage_df = zero_df.merge(coverage_df, left_index=True, right_index=True, how='outer')
    # depth_x "0" column no longer needed.
    coverage_df = coverage_df.drop(columns=['depth_x'])
    coverage_df = coverage_df.rename(columns={'depth_y': 'depth'})
    # Covert the NaN to 0 coverage and get some metrics.
    coverage_df = coverage_df.fillna(0)
    coverage_df['depth'] = coverage_df['depth'].apply(int)
    total_length = len(coverage_df)
    average_coverage = coverage_df['depth'].mean()
    zero_df = coverage_df[coverage_df['depth'] == 0]
    total_zero_coverage = len(zero_df)
    total_coverage = total_length - total_zero_coverage
    genome_coverage = "{:.2%}".format(total_coverage / total_length)
    # Output a zero-coverage vcf fil and the metrics file.
    output_files(vcf_file, total_zero_coverage, zero_df, output_vcf, average_coverage, genome_coverage, output_metrics)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--bam_input', action='store', dest='bam_input', help='bam input file')
    parser.add_argument('--output_metrics', action='store', dest='output_metrics', required=False, default=None, help='Output metrics text file')
    parser.add_argument('--output_vcf', action='store', dest='output_vcf', required=False, default=None, help='Output VCF file')
    parser.add_argument('--reference', action='store', dest='reference', help='Reference dataset')
    parser.add_argument('--vcf_input', action='store', dest='vcf_input', help='vcf input file')

    args = parser.parse_args()

    get_coverage_and_snp_count(args.bam_input, args.vcf_input, args.reference, args.output_metrics, args.output_vcf)
