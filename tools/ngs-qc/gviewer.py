#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on 28/01/15

@author: Matthias Blum
"""

import argparse
import configparser
import libngs
import libregion
import os
import re
import reporting
import shutil
import sys
import tempfile

config = configparser.ConfigParser()
config.read(os.path.join(os.path.dirname(__file__), 'config.ini'))
gene_db = os.path.join(os.path.dirname(__file__), 'genes.db')


#todo: add exiterrors to galaxy_gviewer.xml


def gviewer(infile, loqctable, dest, genome, **kwargs):
    galaxy = kwargs.get('galaxy')
    genes = kwargs.get('genes', [])
    issorted = kwargs.get('issorted', False)
    nodup = kwargs.get('nodup', False)
    nregions = kwargs.get('nregions', 10)
    quiet = kwargs.get('quiet', False)
    regions = kwargs.get('regions', [])
    resolution = kwargs.get('resolution', 500000)
    tmpdir = kwargs.get('tmp')

    if tmpdir is None or not os.path.isdir(tmpdir):
        tmpdir = tempfile.mkdtemp()
        rmtmp = True
    else:
        rmtmp = False

    files = [infile]
    exitcode = 0        # Default exitcode, used if everything goes fine
    libngs.echo(quiet, 'Input file: {}'.format(os.path.basename(infile)))

    bedtools = config.get('softwares', 'bedtools')
    bedutils = os.path.join(os.path.dirname(__file__), 'utils', 'NGSQC_bed_utils')
    gnuplot = config.get('softwares', 'gnuplot')
    wigit = os.path.join(os.path.dirname(__file__), 'utils', 'wigit')
    sort = config.get('softwares', 'sort')
    ncores = config.getint('performance', 'cores')
    maxmem = config.getint('performance', 'memory')

    # Take the input file and get a sorted BED file ready
    # If obj is a tuple, it's fine. If it's an int, something went wrong
    obj = libngs.prepare_bed(infile, galaxy, tmpdir, quiet=quiet, bedtools=bedtools, bedutils=bedutils,
                             issorted=issorted, sort=sort, ncores=ncores, maxmem=maxmem, nodup=nodup)

    if isinstance(obj, tuple):
        # Unpack obj tuple
        files, nreads, nureads, chrms, len_mean = obj

        if regions:
            # Genomic region submitted
            regions2 = []
            for r in regions:
                if r[0] not in chrms:
                    if r[0][3:] not in chrms:
                        r[0] = r[0][3:]
                    else:
                        continue

                regions2.append({
                    'chrm': r[0],
                    'start': r[1],
                    'end': r[2],
                    'genes': libregion.find_genes(genome, r[0], r[1], r[2])
                })

            if regions2:
                regions = regions2
            else:
                sys.stderr.write('Error: chromosome naming problem.\n')
                exitcode = 9

        elif genes:
            # Genes submitted
            regions = []
            for g in genes:
                r = libregion.make_region_around_gene(genome, g, resolution)

                if r:
                    # Verify if the region's chromosome was found in the BED file
                    if r['chrm'] not in chrms:
                        # Dammit!
                        if r['chrm'][3:] in chrms:
                            # Transform "chr10" to "10"
                            r['chrm'] = r['chrm'][3:]
                        else:
                            # Don't support other cases
                            continue

                    regions.append(r)

            if not regions:
                sys.stderr.write('Error: none of the submitted genes were found\n'
                                 'Please, do not use aliases.\n')
                exitcode = 10
        else:
            # Get the best regions

            # Get the size of each chromosome and the (temporary) genome file
            genomeinfo, genomefile = libngs.get_genomeinfo(genome, chrms=chrms)

            if genomeinfo is not None:
                regions, tot_reads = libregion.select_regions(loqctable, genomeinfo, resolution=resolution)

                if regions:
                    # Filter by reads (standard deviation)
                    # Returns a dictionary of regions with the chromosome as key
                    regions = libregion.sortregions(regions, tot_reads)

                    for chrm, reg_lst in regions.items():
                        regions[chrm] = reg_lst[:nregions*5]

                    regions = libregion.filter_by_overlap(regions)

                    for chrm, reg_lst in regions.items():
                        regions[chrm] = reg_lst[:nregions]

                    regions = libregion.filter_by_gene(regions, genome, resolution=resolution,
                                                       minsize=100, mingenes=1)

                    regions = libregion.enqueue_regions(regions, nregions)
                    regions = regions[:nregions]
            else:
                genomeinfo, genomefile = libngs.get_genomeinfo(genome)
                error_msg = ("Error: could not find any reference chromosome in the input BED file.\n"
                             "Please, check the chromosomes "
                             "are properly written (e.g. {} for {}).\n".format(sorted(genomeinfo.keys())[0], genome))
                sys.stderr.write(error_msg)
                exitcode = 8

        if exitcode == 0:
            if regions:
                regions = libregion.get_disp(loqctable, regions)

                regions = libregion.get_wigs(wigit, regions, files[-1])

                libregion.colorbar(os.path.join(tmpdir, 'colorbar.png'), gnuplot=gnuplot)

                imgs = []
                for i, r in enumerate(regions):
                    png = os.path.join(tmpdir, 'region{}.png'.format(i+1))
                    # Show detailed genes if the user submitted a list of genes
                    # or if the resolution is lower or equal to 500kb
                    details = True if genes else r['end'] - r['start'] <= 500000

                    c = libregion.plot(r, png, genedetails=details, gnuplot=gnuplot)
                    if c != 0:
                        break
                    imgs.append(png)

                reporting.gv_report(dest, tmpdir)

                for f in imgs:
                    if os.path.isfile(f):
                        os.unlink(f)
                os.unlink(os.path.join(tmpdir, 'colorbar.png'))
            else:
                sys.stderr.write('Warning: could not select any genomic region\n')
                exitcode = 11
    else:
        # Is an int: the error code. Something went wrong
        exitcode = obj

    libngs.cleanfiles(files)
    if rmtmp:
        shutil.rmtree(tmpdir)

    return exitcode


def check_regions(regions):
    ok_regions = []
    for string in regions:
        m = re.match(r'(.+):(\d+)\-(\d+$)', string, re.I)
        if m:
            chrm = m.group(1)
            start = int(m.group(2))
            end = int(m.group(3))

            if start < end and 10000 <= end - start <= 2000000:
                ok_regions.append([chrm, start, end])

    if not ok_regions:
        sys.stderr.write('Error: invalid submitted regions\n')
        exit(-1)

    return ok_regions


def main():
    parser = argparse.ArgumentParser(description="Display the read count intensity, the bins' dispersion "
                                                 "and the genes of one or more genomic region(s).")

    # Required arguments
    parser.add_argument('input', help='Input BED/BAM file')
    parser.add_argument('localqc', help='Local QC file (10%%)')
    parser.add_argument('-g', dest='genome', required=True,
                        choices=sorted(os.listdir(os.path.join(os.path.dirname(__file__), 'genomes'))),
                        help='Genome assembly')
    parser.add_argument('-o', dest='output', help='PDF output report', required=True)

    # Optional arguments
    parser.add_argument('--genes', nargs='*', help='List of genes')
    parser.add_argument('--nodup', action='store_true', default=False, help='Remove PCR duplicate reads')
    parser.add_argument('--nregions', type=int, default=10, help='Number of local QC regions (default: 10)')
    parser.add_argument('--quiet', action='store_true', default=False, help='Do not print information message')
    parser.add_argument('--regions', nargs='*', help='List of genomic regions (<chr>:<start>-<end>)')
    parser.add_argument('--resolution', type=int, default=500000, help='Resolution (default: 500,000 bp)')
    parser.add_argument('--sorted', action='store_true', default=False, help='Consider the input file sorted')
    parser.add_argument('--tmp', help='Temporary directory')

    args = parser.parse_args()

    # Verify input file's format
    if not args.input.lower().endswith(('.bed.gz', '.bed', '.bam')):
        sys.stderr.write('Error: input file format not supported (only BED/BED.GZ/BAM)\n')
        exit(-1)

    if args.nregions < 0:
        sys.stderr.write('Invalid number of regions ({})\n'.format(args.num_regions))
        exit(-1)

    if args.resolution < 0:
        sys.stderr.write('Invalid resolution ({0})\n'.format(args.resolution))
        exit(-1)

    if args.regions:
        regions = check_regions(args.regions)
    else:
        regions = []

    gviewer(args.input, args.localqc, args.output, args.genome, genes=args.genes,
            issorted=args.sorted, nodup=args.nodup, nregions=args.nregions, quiet=args.quiet, regions=regions,
            resolution=args.resolution, tmp=args.tmp)


if __name__ == '__main__':
    main()


