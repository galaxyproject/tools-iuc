#!/usr/bin/env python

import argparse
import os

from nanocompore.SampCompDB import SampCompDB


def is_valid_file(file_name):
    if os.path.isfile(file_name):
        return os.path.abspath(file_name)
    else:
        raise FileNotFoundError(os.path.abspath(file_name))


def is_valid_directory(dir_name):
    if os.path.isdir(dir_name):
        return os.path.abspath(dir_name)
    else:
        raise NotADirectoryError(os.path.abspath(dir_name))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='save nanocompre sampcomp \
            results as interval outputs \
            \nSample call: \"python Nannocompore-plot.py --db-path \
            ./out_SampComp.db --ref-fasta ref.fa --annotation-bed annot.bed \
            --out-dir ./plots/')

    parser.add_argument('--ref-fasta', required=True, type=is_valid_file,
                        help='The reference genome  used for read alignment.')
    parser.add_argument('--db-path', default="./out_SampComp.db", type=str,
                        help='Path to the SampCompDB database path prefix.')
    parser.add_argument('--annotation-bed', required=False, type=is_valid_file,
                        help='BED file containing the annotation of the transcriptome used as reference when mapping')
    parser.add_argument('--pvalue-types', type=str,
                        default='GMM_logit_pvalue,KS_dwell_pvalue,KS_intensity_pvalue',
                        help='path to the annotations')
    parser.add_argument('--bedgraph', default=False,
                        help='write output in BEDGRAPH format instead of BED')
    parser.add_argument('--pvalue-threshold', default=1.0,
                        help='Maximum reported p-value.')
    parser.add_argument('--out-dir', default="./", type=is_valid_directory,
                        help='path the plotting output directory.')


    args = parser.parse_args()

    db = SampCompDB(args.db_path, fasta_fn=args.ref_fasta,
                    bed_fn=args.annotation_bed)
    print(db)
    print("DB read ids:", db.ref_id_list)

    if args.annotation_bed:
        for pt in args.pvalue_types.split(','):
            print("bedgraph output for p-value type:", pt)
            db.save_to_bed(output_fn='{}/{}.bedgraph'.format(args.out_dir, pt),
                           pvalue_field=pt, pvalue_thr=args.pvalue_threshold,
                           bedgraph=args.bedgraph)
