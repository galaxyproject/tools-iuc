#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on 26/03/15

@author: Matthias Blum
"""

import argparse
import configparser
import libngs
import libquartiles
import libregion
import libtargets
import os
import reporting
import shutil
import sys
import tempfile
import zipfile

from datetime import datetime
from subprocess import Popen, PIPE

VERSION = '1.7.6'

config = configparser.ConfigParser()
config.read(os.path.join(os.path.dirname(__file__), 'config.ini'))


def compute_indicators(bedutils, bed, outdir, genomefile, bg=1, rep=1, gnuplot='gnuplot', binsize=500):
    """
    Divide the input BED file into fixed-size bins and attribute each read to a bin,
    write the supplementary files, then compute the localQC indicators
    :param bedutils:    bedutils binary
    :param bed:         BED file
    :param outdir:      Output directory
    :param genomefile:  Genome file
    :param bg:          Background threshold
    :param rep:         Replicate number
    :param gnuplot:     Gnuplot binary
    :param binsize:     Size of each bin
    :return:
    """
    bins_table = os.path.join(outdir, 'table_samples_replicate_{}.txt.gz'.format(rep))

    # Create a sub-directory for wiggle files
    wigsdir = os.path.join(outdir, 'wigs')
    if not os.path.isdir(wigsdir):
        os.mkdir(wigsdir)

    localqc_table = os.path.join(outdir, 'sam_10pc_all_replicate_{}.txt'.format(rep))
    bed_track = os.path.join(outdir, 'local_QC-scores_replicate_{}.s50_10pc.bed'.format(rep))
    indicators_file = os.path.join(outdir, 'indicators.txt')

    cmd = [bedutils, 'bin_table', bed, genomefile, str(binsize), str(bg), outdir, '--plot',
           '--gnuplot-path', gnuplot, '--globalQC', '--sam-table', localqc_table,
           '--bed-track', bed_track, '--wigs', wigsdir]

    with open(bins_table, 'wb') as fh, open(os.devnull, 'w') as dn:
        pop1 = Popen(cmd, stdout=PIPE, stderr=dn)
        pop2 = Popen(['gzip', '-c'], stdin=pop1.stdout, stdout=fh, stderr=dn)
        pop2.wait()
        pop1.wait()

    files = [localqc_table, bed_track, bins_table, indicators_file]
    efiles = [f for f in files if os.path.isfile(f)]

    if len(files) != len(efiles) or pop2.returncode != 0:
        print("{} {} {}".format(files,efiles,pop2.returncode))
        for f in efiles:
            os.unlink(f)
        return None, None

    indicators = {}
    with open(indicators_file) as f:
        for line in f:
            k, v = line.rstrip().split('\t', 1)
            if k == 'Bins':
                indicators['bins'] = int(v)
            elif k == 'Bins w/o background':
                indicators['bins2'] = int(v)
            elif k == 'DenQC_s90_2.5pc':
                indicators['denqc_90_2.5'] = float(v)
            elif k == 'DenQC_s90_5pc':
                indicators['denqc_90_5'] = float(v)
            elif k == 'DenQC_s90_10pc':
                indicators['denqc_90_10'] = float(v)
            elif k == 'DenQC_s50_2.5pc':
                indicators['denqc_50_2.5'] = float(v)
            elif k == 'DenQC_s50_5pc':
                indicators['denqc_50_5'] = float(v)
            elif k == 'DenQC_s50_10pc':
                indicators['denqc_50_10'] = float(v)
            elif k == 'SimQC_2.5pc':
                indicators['simqc_2.5'] = float(v)
            elif k == 'SimQC_5pc':
                indicators['simqc_5'] = float(v)
            elif k == 'SimQC_10pc':
                indicators['simqc_10'] = float(v)

    # We don't need this file any more
    os.unlink(indicators_file)

    # Renaming files
    os.rename(os.path.join(outdir, 'pc_s50_s70_s90_bgs.png'),
              os.path.join(outdir, 'pc_s50_s70_s90_replicate_{}.png'.format(rep)))

    os.rename(os.path.join(outdir, 'dispersion_s50_s70_s90_bgs.png'),
              os.path.join(outdir, 'dispersion_s50_s70_s90_replicate_{}.png'.format(rep)))

    for f in os.listdir(wigsdir):
        if f[-4:] == '.wig':
            chrm, string = f.split('_', 1)
            os.rename(os.path.join(wigsdir, f),
                      os.path.join(wigsdir, '{}_replicate_{}.{}'.format(chrm, rep, string)))

    return bins_table, indicators


def plot_regions(bed, localqc_table, genome, genomeinfo, rep, outdir, wigit, gnuplot, bg=1,
                 nregions=10, lowres=500000, highres=50000):
    """
    Select and plot a given number of genomic regions, based on the localQC scores,
    the read count intensity, or the presence of genes
    :param bed:             Input BED file
    :param localqc_table:   LocalQC table
    :param genome:          Genome assembly
    :param genomeinfo:      Dictionary of chromosome size
    :param rep:             Current replicate
    :param outdir:          Output directory
    :param wigit:           Wigit binary
    :param gnuplot:         Gnuplot binary
    :param bg:              Background threshold
    :param nregions:        Number of regions to select
    :param lowres:          Low resolution
    :param highres:         High resolution
    :return:                True if regions were selected/plotted, False if we did not find any region or if
    """
    regions, tot_reads = libregion.select_regions(localqc_table, genomeinfo, resolution=lowres)

    if regions:
        # Filter by reads (standard deviation)
        # Returns a dictionary of regions with the chromosome as key
        regions = libregion.sortregions(regions, tot_reads)

        # Since we have a huuuuge list of regions, discard some
        for chrm, reg_lst in regions.items():
            regions[chrm] = reg_lst[:nregions*5]

        # Remove overlapping regions (more than 25%)
        regions = libregion.filter_by_overlap(regions)

        # Keep the total number of regions we want, but for EACH chromosomes
        # That should be enough to have nice regions without slowing down SQlite
        for chrm, reg_lst in regions.items():
            regions[chrm] = reg_lst[:nregions]

        # Keep regions with at least one gene
        regions = libregion.filter_by_gene(regions, genome, resolution=highres,
                                           minsize=100, mingenes=1)

        # Keep the best region of each chromosome,
        # and the best regions (all chromosome confounded)
        regions = libregion.enqueue_regions(regions, nregions)

        regions = regions[:nregions]

        # Get the dispersion of each local QC bin in each region
        regions = libregion.get_disp(localqc_table, regions, bg=bg)

        # Get the wiggles for each region
        regions = libregion.get_wigs(wigit, regions, bed)

        libregion.colorbar(os.path.join(outdir, 'colorbar.png'), gnuplot=gnuplot)

        for i, r in enumerate(regions):
            png = os.path.join(outdir, 'localqc.rep{}.region{}.png'.format(rep, i * 2 + 1))
            c = libregion.plot(r, png, genedetails=False, gnuplot=gnuplot)

            if c != 0:
                return False

            r = libregion.shrink_region(r, resolution=highres)

            png = os.path.join(outdir, 'localqc.rep{}.region{}.png'.format(rep, i * 2 + 2))
            c = libregion.plot(r, png, genedetails=True, gnuplot=gnuplot)

            if c != 0:
                return False
    else:
        sys.stderr.write('Warning: could not select any genomic region\n')
        return False

    return True


def write_summary(dest, filename, nreads, nureads, mean_len, rep, nrep, indicators, chrms, genome,
                  target=None, nodup=False, nobgs=True, binsize=500, bg=1, ci=0.995):
    """
    Write a textual summary file for the processed run
    :param dest:        Output file
    :param filename:    Input file name
    :param nreads:      Number of total reads
    :param nureads:     Number of unique reads
    :param mean_len:    Mean of reads' size
    :param rep:         Current replicate
    :param nrep:        Number of replicates
    :param indicators:  QC indicators
    :param chrms:       List of chromosomes found in the BED file
    :param genome:      Genome assembly
    :param target:      Target molecule
    :param nodup:       PCR dupe we removed
    :param nobgs:       Background wasn't subtracted
    :param binsize:     Size of each bin
    :param bg:          Background threshold
    :param ci:          Confidence interval
    :return:
    """
    with open(dest, 'w') as fh:
        fh.write('Generated by NGS-QC Generator {0} on {1}\n'.format(VERSION,
                                                                     datetime.now().strftime('%Y-%m-%d %H:%M')))
        fh.write('Clonal reads removal\t{}\n'.format('On' if nodup else 'Off'))
        fh.write('Background subtraction\t{}\n'.format('Off' if nobgs else 'On'))
        fh.write('Input file\t{}\n'.format(filename))
        fh.write('Genome assembly\t{}\n'.format(genome))
        fh.write('Target molecule\t{}\n'.format(target))
        fh.write('Replicate\t{}/{}\n'.format(rep, nrep))
        fh.write('Windows size\t{}\n'.format(binsize))
        fh.write('Sampling\t{}\n'.format(','.join(str(i) for i in (0.5, 0.7, 0.9))))
        fh.write('Reads\t{}\n'.format(nreads))
        fh.write('Unique reads\t{}\n'.format(nureads))
        fh.write("Reads size mean\t{:.2f}\n".format(mean_len))
        fh.write('Bins\t{}\n'.format(indicators['bins']))
        fh.write('DenQC_s90_2.5pc(%)\t{}\n'.format(indicators['denqc_90_2.5']))
        fh.write('DenQC_s50_2.5pc(%)\t{}\n'.format(indicators['denqc_50_2.5']))
        fh.write('SimQC_2.5pc (DenQC_s90/s50)\t{}\n'.format(indicators['simqc_2.5']))
        fh.write('DenQC_s90_5pc(%)\t{}\n'.format(indicators['denqc_90_5']))
        fh.write('DenQC_s50_5pc(%)\t{}\n'.format(indicators['denqc_50_5']))
        fh.write('SimQC_5pc (DenQC_s90/s50)\t{}\n'.format(indicators['simqc_5']))
        fh.write('DenQC_s90_10pc(%)\t{}\n'.format(indicators['denqc_90_10']))
        fh.write('DenQC_s50_10pc(%)\t{}\n'.format(indicators['denqc_50_10']))
        fh.write('SimQC_10pc (DenQC_s90/s50)\t{}\n'.format(indicators['simqc_10']))
        fh.write('Chromosomes\t{}\n'.format(','.join(chrms)))
        fh.write('Threshold\t{}\n'.format(bg))
        fh.write('Bins w/o background\t{}\n'.format(indicators['bins2']))
        fh.write('Conf. interval\t{}\n'.format(ci))


def ngsqc(infile, outdir, genome, **kwargs):
    """
    NGS-QC Generator function: process a input file
    :param infile:  Input BED[.GZ]/BAM file
    :param outdir:  Results output directory
    :param genome:  Genome assembly
    :param kwargs:
    :return:
    """
    galaxy = kwargs.get('galaxy')               # Input file name for Galaxy
    highres = kwargs.get('highres', 50000)      # High resolution
    issorted = kwargs.get('issorted', False)    # Consider the input file sorted
    lowres = kwargs.get('lowres', 500000)       # Low resolution for genomic regions
    nobgs = kwargs.get('nobgs', False)          # Do not subtract the background
    nodup = kwargs.get('nodup', False)          # Remove the PCR duplicate reads
    nregions = kwargs.get('regions', 10)        # Number of genomic regions to display
    nreplicates = kwargs.get('replicates', 1)   # Number of virtual replicates
    quiet = kwargs.get('quiet', False)          # Do not show information messages (errors are still showed)
    target = kwargs.get('target')               # Target molecule
    tmpdir = kwargs.get('tmp')                  # Temporary directory

    if tmpdir is None:
        tmpdir = tempfile.mkdtemp()
    else:
        tmpdir = os.path.join(tmpdir, 'ngsqc' + os.path.basename(infile))

        if not os.path.isdir(tmpdir):
            os.makedirs(tmpdir)

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    files = [infile]
    exitcode = 0        # Default exitcode, used if everything goes fine

    libngs.echo(quiet, 'Input file: {}'.format(os.path.basename(infile)))

    bedtools = config.get('softwares', 'bedtools')
    bedutils = os.path.join(os.path.dirname(__file__), 'utils', 'NGSQC_bed_utils')
    gnuplot = config.get('softwares', 'gnuplot')
    wigit = os.path.join(os.path.dirname(__file__), 'utils', 'wigit', 'wigit')
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

        # Convert all paths to real paths
        files = [os.path.realpath(f) for f in files]

        # Get the size of each chromosome and the (temporary) genome file
        genomeinfo, genomefile = libngs.get_genomeinfo(genome, chrms=chrms)

        if genomeinfo is not None:
            if nobgs:
                bg = 1
            elif nodup:
                # Considered reads: unique reads if we removed the PCR dupe reads
                bg = libngs.getbackground(nureads, sum(genomeinfo.values()), binsize=500, ci=0.995)
            else:
                bg = libngs.getbackground(nreads, sum(genomeinfo.values()), binsize=500, ci=0.995)

            if bg is not None:
                for rep in range(1, nreplicates + 1):
                    libngs.echo(quiet, '\tReplicate #{}: computing indicators'.format(rep))

                    # Compressed bin table, directory with the supplementary files, indicators dictionary
                    bins_table, indicators = compute_indicators(bedutils, files[-1], tmpdir, genomefile,
                                                                bg=bg, rep=rep, gnuplot=gnuplot, binsize=500)

                    if bins_table is not None:
                        localqcdir = os.path.join(tmpdir, 'localqc')
                        if nregions > 0:
                            if not os.path.isdir(localqcdir):
                                os.mkdir(localqcdir)

                            # Give the proper order, even if the passed values are messed up (high > low)
                            highres, lowres = sorted([lowres, highres])
                            localqc_table = os.path.join(tmpdir, 'sam_10pc_all_replicate_{}.txt'.format(rep))

                            plot_regions(files[-1], localqc_table, genome, genomeinfo, rep, localqcdir,
                                         bg=bg, wigit=wigit, gnuplot=gnuplot,
                                         nregions=nregions, lowres=lowres, highres=highres)

                        filename = os.path.basename(infile) if galaxy is None else galaxy
                        summary = os.path.join(tmpdir, 'summary_replicate_{}.txt'.format(rep))
                        write_summary(summary, filename, nreads, nureads, len_mean, rep, nreplicates, indicators,
                                      chrms, genome, target=target, nodup=nodup,
                                      nobgs=nobgs, binsize=500, bg=bg, ci=0.995)

                        dbquartiles, dbversion = libquartiles.get()
                        if target:
                            organism = libngs.getorganism(genome)
                            public_data, target = libtargets.get_data(target, organism)
                            organism_abbr = ''.join(i[0].lower() for i in organism.split())
                            public_label = 'public {} ({})'.format(target, organism_abbr)

                            qcvals = {
                                2.5: libngs.getqcval(indicators['denqc_50_2.5'], indicators['simqc_2.5']),
                                5: libngs.getqcval(indicators['denqc_50_5'], indicators['simqc_5']),
                                10: libngs.getqcval(indicators['denqc_50_10'], indicators['simqc_10'])
                            }

                            for disp in (2.5, 5, 10):
                                dest = os.path.join(tmpdir, 'public.{}.png'.format(disp))
                                # public_data contains qcvalues for all dispersion:
                                # we only want the qcvalue for the current disp
                                data_disp = [(p_reads, p_qcvals[disp]) for p_reads, p_qcvals in public_data]
                                libtargets.plot(data_disp, disp, nreads, qcvals[disp],
                                                dbquartiles[disp], dest, label=target, public_label=public_label)

                        pdf = os.path.join(tmpdir, 'NGS-QC_report_replicate_{}.pdf'.format(rep))
                        r = reporting.init(summary, pdf, tmpdir, quartiles=dbquartiles, dbversion=dbversion)

                        if r is None:
                            sys.stderr.write('Error while creating report\n')
                            exitcode = 11
                        else:
                            # Move PDF report and bin table
                            if os.path.realpath(tmpdir) != os.path.realpath(outdir):
                                pdf_dest = os.path.join(outdir, os.path.basename(pdf))
                                if os.path.isfile(pdf_dest):
                                    os.unlink(pdf_dest)
                                shutil.move(pdf, pdf_dest)

                                table_dest = os.path.join(outdir, os.path.basename(bins_table))
                                if os.path.isfile(table_dest):
                                    os.unlink(table_dest)
                                shutil.move(bins_table, table_dest)

                            # Compress supplementary files
                            zippath = os.path.join(outdir, 'local_QC_indicators_replicate_{}.zip'.format(rep))
                            with zipfile.ZipFile(zippath, 'w', compression=zipfile.ZIP_DEFLATED) as z:
                                for f in os.listdir(tmpdir):
                                    path = os.path.realpath(os.path.join(tmpdir, f))

                                    # The temporary BED file is in tmpdir, but it must not be included in the zip file!
                                    if os.path.isfile(path) and path not in files:
                                        z.write(path, f)
                                        os.unlink(path)

                                for f in os.listdir(localqcdir):
                                    path = os.path.join(localqcdir, f)
                                    z.write(path, os.path.join('localqc', f))
                                    os.unlink(path)

                                for f in os.listdir(os.path.join(tmpdir, 'wigs')):
                                    path = os.path.join(os.path.join(tmpdir, 'wigs'), f)
                                    if os.path.isfile(path):
                                        z.write(path, os.path.join('wigs', f))
                                        os.unlink(path)

                        if os.path.isdir(localqcdir):
                            shutil.rmtree(localqcdir)

                        if os.path.isdir(os.path.join(tmpdir, 'wigs')):
                            shutil.rmtree(os.path.join(tmpdir, 'wigs'))

                        if exitcode:
                            break
                    else:
                        error_msg = ('Error while computing indicators.\n'
                                     'You might have selected the wrong genome assembly or '
                                     'one or more reads are outside of the reference genome.\n')
                        sys.stderr.write(error_msg)

                        # If a replicate fails, we don't process the others
                        exitcode = 10
                        break
            else:
                sys.stderr.write('Error during the Poisson distribution computation\n')
                exitcode = 9

            os.unlink(genomefile)
        else:
            genomeinfo, genomefile = libngs.get_genomeinfo(genome)
            error_msg = ("Error: could not find any reference chromosome in the input BED file.\n"
                         "Please, check the chromosomes "
                         "are properly written (e.g. {} for {}).\n".format(sorted(genomeinfo.keys())[0], genome))
            sys.stderr.write(error_msg)
            exitcode = 8
    else:
        # Is an int: the error code. Something went wrong
        exitcode = obj

    libngs.cleanfiles(files)
    shutil.rmtree(tmpdir)

    return exitcode


def main():
    parser = argparse.ArgumentParser(description='NGS-QC Generator: Control the quality of profiles'
                                                 'obtained by ChIP-sequencing')

    # Required arguments
    parser.add_argument('input', help='Input BED/BAM file')
    parser.add_argument('-o', dest='outdir', default=os.getcwd(), help='Output directory (default: working directory)')
    parser.add_argument('-g', dest='genome', required=True,
                        choices=sorted(os.listdir(os.path.join(os.path.dirname(__file__), 'genomes'))),
                        help='Genome assembly')

    # Optional arguments
    parser.add_argument('--highres', type=int, default=50000, help='LocalQC: high resolution (default: 50,000 bp)')
    parser.add_argument('--lowres', type=int, default=500000, help='LocalQC: low resolution (default: 500,000 bp)')
    parser.add_argument('--nobgs', action='store_true', default=False, help='Do not remove the background')
    parser.add_argument('--nodup',  action='store_true', default=False, help='Remove PCR duplicate reads')
    parser.add_argument('--quiet', action='store_true', default=False, help='Do not print information message')
    parser.add_argument('--replicates', type=int, default=1, help='Number of virtual replicates (default: 1)')
    parser.add_argument('--regions', type=int, default=10, help='Number of local QC regions (default: 10)')
    parser.add_argument('--sorted', action='store_true', default=False, help='Consider the input file sorted')
    parser.add_argument('--target', help='Input target molecule')
    parser.add_argument('--tmp', help='Temporary directory')

    args = parser.parse_args()

    # Verify input file's format
    if not args.input.lower().endswith(('.bed.gz', '.bed', '.bam')):
        sys.stderr.write('Error: input file format not supported (only BED/BED.GZ/BAM)\n')
        exit(-1)

    if args.regions < 0:
        args.regions = 0

    if args.lowres < 100000:
        args.lowres = 100000

    if args.highres < 10000:
        args.highres = 10000

    if args.replicates < 0:
        sys.stderr.write('Error: invalid number of replicates ({})\n'.format(args.replicates))
        exit(-1)

    if args.tmp is not None and not os.path.isdir(args.tmp):
        sys.stderr.write('Error: temporary directory "{}" does not exist\n'.format(args.tmp))
        exit(-1)

    ngsqc(args.input, args.outdir, args.genome, regions=args.regions, lowres=args.lowres,
          highres=args.highres, replicates=args.replicates, tmp=args.tmp, nobgs=args.nobgs,
          nodup=args.nodup, issorted=args.sorted,
          quiet=args.quiet, target=args.target)


if __name__ == "__main__":
    main()
