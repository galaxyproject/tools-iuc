#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on 27/03/15

@author: Matthias Blum
"""

import argparse
import os
import shutil
import time
import zipfile
from datetime import datetime

import gviewer
import ngsqc


def write_html(filename, genome, files, dest, exec_time):
    """
    Write an HTML page for Galaxy
    :param filename:    Input file naùe
    :param genome:      Genome assembly
    :param files:       Output files (report + zip)
    :param dest:        HTML output file
    :param exec_time:   Execution time
    :return:
    """
    html_string = ('<html><head><style>'
                   'table {{border-collapse: collapse;margin: 0 auto;}}'
                   'td {{padding: 5px;border: 1px solid #000;border-collapse: collapse;}}'
                   'th.a {{width: 200px;}}'
                   'h2 {{margin: 0;}}'
                   'td.b {{text-align: center;background-color: #396A92; color: #FFFFFF;}}'
                   '</style></head>'
                   '<body><table><thead><tr><th class="a"></th><th colspan="2"></th></tr></thead>'
                   '<tbody>'
                   '<tr><td colspan="3" class="b">'
                   '<h2>NGS-QC Generator results</h1></td>'
                   '<tr><td>Input file</td>'
                   '<td colspan="2">{}</td></tr>'
                   '<tr><td>Genome assembly</td>'
                   '<td colspan="2">{}</td></tr>'
                   '<tr><td>Page generated on</td>'
                   '<td colspan="2">{}</td></tr>'
                   "<tr><td>Run time</td>"
                   "<td colspan='2'>{:.2f} seconds</td></tr>").format(filename,
                                                                      genome,
                                                                      datetime.now().strftime('%Y-%m-%d %H:%H'),
                                                                      exec_time)

    for i, rep_files in enumerate(files):
        pdf = os.path.basename(rep_files[0])
        zipf = os.path.basename(rep_files[1])
        html_string += ("<tr><td>Replicate {}</td>"
                        "<td><a href=\"{}\">PDF report</a></td>"
                        "<td><a href=\"{}\">Supplementary data</a></td></tr>".format(i + 1, pdf, zipf))

    html_string += "</tbody></table></body></html>"

    with open(dest, 'w') as fh:
        fh.write(html_string)


def main():
    parser = argparse.ArgumentParser(description='Galaxy wrapper for NGS-QC tools')

    subparsers = parser.add_subparsers(title='command', dest='command')

    # NGS-QC
    parser1 = subparsers.add_parser('ngsqc', help='NGS-QC Generator')
    parser1.add_argument('input', help='Input BED/BAM file')
    parser1.add_argument('-o', dest='outdir', help='Output directory')
    parser1.add_argument('-g', dest='genome', required=True,
                         choices=sorted(os.listdir(os.path.join(os.path.dirname(__file__), 'genomes'))),
                         help='Genome assembly')
    parser1.add_argument('--nobgs', action='store_true', default=False, help='Do not remove the background')
    parser1.add_argument('--nodup',  action='store_true', default=False, help='Remove PCR duplicate reads')
    parser1.add_argument('--lowres', type=int, default=500000, help='LocalQC: low resolution (default: 500,000 bp)')
    parser1.add_argument('--replicates', type=int, default=1, help='Number of virtual replicates (default: 1)')
    parser1.add_argument('--target', help='Input target molecule')
    parser1.add_argument('--filename', help='Galaxy dataset')
    parser1.add_argument('--tables', nargs='+', help='Output localQC (10%%) tables (galaxy)')
    parser1.add_argument('--html', help='Output HTML file (galaxy)')

    # Genome Viewer
    parser2 = subparsers.add_parser('gviewer', help='Genome viewer')
    parser2.add_argument('input', help='Input BED/BAM file')
    parser2.add_argument('localqc', help='Local QC file (10%%)')
    parser2.add_argument('-g', dest='genome', required=True,
                         choices=sorted(os.listdir(os.path.join(os.path.dirname(__file__), 'genomes'))),
                         help='Genome assembly')
    parser2.add_argument('-o', dest='output', help='PDF output report', required=True)
    parser2.add_argument('--filename', help='Galaxy dataset')
    parser2.add_argument('--genes', nargs='*', help='List of genes')
    parser2.add_argument('--nodup', action='store_true', default=False, help='Remove PCR duplicate reads')
    parser2.add_argument('--nregions', type=int, default=10, help='Number of local QC regions (default: 10)')
    parser2.add_argument('--regions', nargs='*', help='List of genomic regions (<chr>:<start>-<end>)')
    parser2.add_argument('--resolution', type=int, default=500000, help='Resolution (default: 500,000 bp)')

    args = parser.parse_args()

    if args.command == 'ngsqc':
        if args.target == 'Unspecified':
            args.target = None

        exec_time = time.time()
        c = ngsqc.ngsqc(args.input, args.outdir, args.genome, nobgs=args.nobgs, lowres=args.lowres,
                        replicates=args.replicates, nodup=args.nodup,
                        target=args.target, quiet=True, galaxy=args.filename)
        exec_time = time.time() - exec_time

        if c == 0:
            files = []
            for rep in range(1, args.replicates + 1):
                pdf = os.path.join(args.outdir, 'NGS-QC_report_replicate_{}.pdf'.format(rep))
                zipf = os.path.join(args.outdir, 'local_QC_indicators_replicate_{}.zip'.format(rep))
                localqc_table = 'sam_10pc_all_replicate_{}.txt'.format(rep)

                files.append([pdf, zipf])

                with zipfile.ZipFile(zipf, 'r', compression=zipfile.ZIP_DEFLATED) as z:
                    z.extract(localqc_table.format(rep), path=args.outdir)

                shutil.move(os.path.join(args.outdir, localqc_table), args.tables[rep - 1])

            write_html(args.filename, args.genome, files=files, dest=args.html, exec_time=exec_time)

        exit(c)
    elif args.command == 'gviewer':
        if args.regions:
            regions = gviewer.check_regions(args.regions[0].split('__cn__'))
        else:
            regions = []

        if args.genes:
            args.genes = args.genes[0].split('__cn__')

        c = gviewer.gviewer(args.input, args.localqc, args.output, args.genome, galaxy=args.filename, genes=args.genes,
                            nodup=args.nodup, nregions=args.nregions, quiet=True, regions=regions,
                            resolution=args.resolution)

        exit(c)


if __name__ == '__main__':
    main()
