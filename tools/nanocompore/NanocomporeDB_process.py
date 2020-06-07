#!/usr/bin/env python

import argparse
import os

from nanocompore.SampCompDB import SampCompDB
# from nanocompore.common import *
# import ast
# from tqdm import tqdm
# from matplotlib import pyplot


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

# def doall(db, reference_fa, gene, rangeL, rangeR, outdir, annot_bed,
#           plot_signals):
#     dict_genes = {}
#     for g in [gene]:
#         lref = [ref for ref in db.ref_id_list if g in ref]
#         assert(len(lref)==1)
#         dict_genes[g]= lref[0]
#     print("genes name map:", dict_genes)
#     for g in dict_genes:
#         fig, ax = db.plot_pvalue (dict_genes[g], tests="logit")
#         fig.savefig(outdir+"/"+g+'-full-pvalue.svg')
#         pyplot.close(fig)
#         fig, ax = db.plot_coverage (dict_genes[g])
#         fig.savefig(outdir+"/"+g+'-full-coverage.svg')
#         pyplot.close(fig)
#     db_entry = db[dict_genes['chrSCV']]
#     for  reg in [[gene, rangeL, rangeR]]:
#         gene, pleft, pright = reg[0], reg[1], reg[2]
#         fig_prefix = outdir+'/{}-{}-{}'.format(gene, pleft, pright)
#         fig, ax = db.plot_coverage (dict_genes[gene], pleft, pright)
#         fig.savefig(fig_prefix + '-coverage.svg')
#         pyplot.close(fig)
#         fig, ax = db.plot_pvalue(dict_genes[gene],pleft, pright,
#                                   tests="logit")
#         fig.savefig(fig_prefix + '-pvalue.svg')
#         pyplot.close(fig)
#         if plot_signals is True:
#             if pright-pleft>200:
#                 pleft, pright = int((pleft+pright)/2-100),
#                                   int((pleft+pright)/2+100)
#                 fig_prefix += '-zoom'
#             fig, ax = db.plot_signal (dict_genes[gene],pleft, pright,
#                                           figsize=(60,10))
#             fig.savefig(fig_prefix +'-signal-violin.svg',
#                           bbox_inches="tight")
#             pyplot.close(fig)
#             if pright - pleft >10:
#                 pleft, pright = int((pleft+pright)/2-2),
#                                   int((pleft+pright)/2 + 2)
# #                 fig_prefix += '-supzoomed'
#             fig, ax = db.plot_signal (dict_genes[gene], pleft, pright,
#                                           kind="swarmplot")
#             fig.savefig(fig_prefix +'-signal-swarm.png', bbox_inches="tight")
#             pyplot.close(fig)


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

    # parser.add_argument('--gene-range', required=True,
    #                   help='The plot gene range, e.g. SOX:1000-2000')
    # parser.add_argument('--plot-signals', default=False,
    #                   help='plot signal violin and swarm plots')

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

    # gene_name = args.gene_range.split(":")[0]
    # rangeL,rangeR = [int(i) for i in args.gene_range.split(":")[1].split("-")]
    # doall(args.db_path, args.ref_fasta, gene_name, rangeL, rangeR,
    #       args.out_dir, args.annotation_bed, args.plot_signals)
