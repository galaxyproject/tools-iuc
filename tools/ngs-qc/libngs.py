#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on 26/03/15

@author: Matthias Blum
"""

import math
import os
import sys
import tempfile
from scipy.stats.distributions import poisson
from subprocess import Popen, PIPE


def getorganism(genome):
    organisms = {
        'Homo sapiens': ['hg18', 'hg19'],
        'Mus musculus': ['mm8', 'mm9', 'mm10'],
        'Caenorhabditis elegans': ['ce4', 'ce6', 'ce10'],
        'Danio rerio': ['danRer5', 'danRer6', 'danRer7'],
        'Drosophila melanogaster': ['dm3'],
        'Rattus norvegicus': ['rn3', 'rn4', 'rn5'],
        'Arabidopsis thaliana': ['TAIR10'],
        'Gallus gallus': ['galGal4']
    }

    for k, v in organisms.items():
        if genome in v:
            return k
    return None


def getqcval(denqc_50, simpqc):
    """
    Compute the QC value for a dispersion's simQC/denQC
    :param denqc_50:
    :param simpqc
    :return:
    """
    try:
        qc_val = math.log(denqc_50/simpqc, 2)
    except ValueError as e:
        qc_val = None
    finally:
        return qc_val


def num_format(x):
    """
    Return a number formated with comma as thousand-separator
    :param x:
    :return:
    """
    return str('{:,}'.format(int(x)))


def bamtobed(path, dest, bedtools='bedtools'):
    with open(dest, 'w') as fh, open(os.devnull, 'w') as dn:
        pop = Popen([bedtools, 'bamtobed', '-i', path], stdout=fh, stderr=dn)
        pop.wait()

    if os.path.isfile(dest):
        if pop.returncode == 0:
            return dest
        os.unlink(dest)
    return None


def cleanfiles(files):
    for f in files:
        if f != files[0] and os.path.isfile(f):
            os.unlink(f)


def echo(quiet, string):
    if not quiet:
        print(string)


def gunzip(path, dest, d=False):
    """

    :param path:
    :param dest:
    :param d:
    :return:
    """
    if d:
        fh = open(dest, 'w')
        cmd = ['gzip', '-cd', path]
    else:
        fh = open(dest, 'wb')
        cmd = ['gzip', '-c', path]

    with open(os.devnull, 'w') as dn:
        pop = Popen(cmd, stdout=fh, stderr=dn)
        pop.wait()

    fh.close()

    if os.path.isfile(dest):
        if pop.returncode == 0:
            return dest
        os.unlink(dest)
    return None


def isbinary(infile):
    """
    Check if a file is a binary file or a text file. Actually, just read the first line.
    :param input_file:
    @return:
    """
    textchars = ''.join(map(chr, [7, 8, 9, 10, 12, 13, 27] + range(0x20, 0x100)))
    isbinary_string = lambda x: bool(x.translate(None, textchars))

    with open(infile) as f:
        b = isbinary_string(f.read(1024))

    return b


def bedintegrity(bedutils, bed):
    pop = Popen([bedutils, 'integrity', bed], stdout=PIPE, stderr=PIPE)
    out, err = pop.communicate()

    if pop.returncode != 0:
        return tuple([None] * 10)

    try:
        cols = out.decode('utf-8').rstrip().split('\t', 9)
        errno = int(cols[0])                    # if equal 0, BED is good
        reads = int(cols[1])                    # Number of reads
        ureads = int(cols[2])                   # Number of considered reads (if a genome file is given)
        reads_length_mean = float(cols[3])      # Mean of reads' length
        dif_read_size = bool(int(cols[4]))      # True if all reads don't have the same length
        min_read_size = int(cols[5])            # Smallest read
        max_read_size = int(cols[6])            # Longest read
        noheader = int(cols[7]) == 0            # True if a no BED header is present
        chrms = sorted(cols[8].split(','))      # List of chromosomes found in the BED file
        last_line = cols[9]                     # Last processed line
    except Exception as e:
        print(e)
        return tuple([None] * 10)
    else:
        return (errno, reads, ureads, reads_length_mean, dif_read_size,
                min_read_size, max_read_size, noheader, chrms, last_line)


def isbedsorted(bedutils, bed, noheader=False):
    """
    Verify if a BED file is sorted
    :param bedutils:    Bedutils binary
    :param bed:         Input BED file
    :param noheader:    True if we don't need to check for BED headers
    :return:            None if an error occured, False if the BED is not sorted, True if it is
    """
    if noheader:
        cmd = [bedutils, 'issorted', bed, '--no-header']
    else:
        cmd = [bedutils, 'issorted', bed]

    pop = Popen(cmd, stdout=PIPE, stderr=PIPE)
    out, err = pop.communicate()

    if pop.returncode != 0:
        return None

    try:
        code = int(out.decode('utf-8').rstrip())
    except Exception as e:
        print("Error to get sorted {}".format(e))
        return None
    else:
        return bool(code)


def count_unique_reads(bedutils, bed, noheader=False):
    if noheader:
        cmd = [bedutils, 'unique', bed, '--no-header']
    else:
        cmd = [bedutils, 'unique', bed]

    pop = Popen(cmd, stdout=PIPE, stderr=PIPE)
    out, err = pop.communicate()

    if pop.returncode != 0:
        return None

    try:
        nunique = int(err.decode('utf-8').rstrip().split('\t')[1])
    except:
        return None
    else:
        return nunique


def filter_unique_reads(bedutils, bed, dest=None, noheader=False):
    if noheader:
        cmd = [bedutils, 'unique', bed, '--write', '--no-header']
    else:
        cmd = [bedutils, 'unique', bed, '--write']

    fh = open(os.devnull, 'w') if dest is None else open(dest, 'w')
    pop = Popen(cmd, stdout=fh, stderr=PIPE)
    out, err = pop.communicate()

    try:
        nunique = int(err.decode('utf-8').rstrip().split('\t')[1])
    except Exception as e:
        nunique = None

    if dest is not None:
        if os.path.isfile(dest) and nunique is not None and pop.returncode == 0:
            return nunique
        elif os.path.isfile(dest):
            os.unlink(dest)
        return None
    elif nunique is not None and pop.returncode == 0:
        return nunique
    else:
        return None


def prepare_bed(infile, galaxy, tmpdir, bedutils, bedtools='bedtools', quiet=False, issorted=False,
                sort='sort', ncores=None, maxmem=None, nodup=False):
    """
    Take an input BED/BAM file and do what it takes to have a proper, sorted BED file
    :param infile:      Input file
    :param galaxy:      Galaxy dataset's filename
    :param tmpdir:      Temporary directory
    :param bedutils:    Bedutils binary
    :param bedtools:    Bedtools binary
    :param quiet:       If True, do not print information message (errors are still printed)
    :param issorted:    Consider the input file sorted
    :param sort:        Unix sort binary
    :param ncores:      Number of cores (sort)
    :param maxmem:      Maximum GB of memory to use (sort)
    :param nodup:       Remove PCR dupe reads
    :return:            An int (error code) if an error occurs, otherwise a tuple
    """
    # Track history of all files we use
    files = [infile]

    # Define the input file's format
    if galaxy:
        # Either BAM or BED. Gzip files are decompressed by Galaxy.
        informat = 'bam' if isbinary(infile) else 'bed'
    elif infile.lower().endswith('.bam'):
        informat = 'bam'
    elif infile.lower().endswith('.bed.gz'):
        informat = 'gz'
    else:
        informat = 'bed'

    if informat == 'bam':
        echo(quiet, '\tBAM to BED conversion')
        bed = os.path.join(tmpdir, os.path.basename(infile)[:-4] + '.bed')
        bed = bamtobed(infile, bed, bedtools)
        if bed is None:
            sys.stderr.write('Error during BAM to BED conversion\n')
            return 1
    elif informat == 'gz':
        echo(quiet, '\tExtracting BED')
        bed = os.path.join(tmpdir, os.path.basename(infile)[:-3])
        bed = gunzip(infile, bed, d=True)
        if bed is None:
            sys.stderr.write('Error during BED file extraction\n')
            return 2
    else:
        bed = infile

    files.append(bed)

    # At this point, we have a BED file
    # Let's check if it's not ill-formated
    echo(quiet, "\tVerifying BED file's integrity")

    errno, nreads, ureads, len_mean, var_len, min_len, max_len, noheader, chrms, last_line = bedintegrity(bedutils, bed)
    if errno is None:
        sys.stderr.write('Error in bedutils (integrity)\n')
        return 3
    elif errno != 0:
        sys.stderr.write('Error: BED file ill-formatted:\n{}\n'.format(last_line))
        return 4

    if not issorted:
        # Verify if the BED file is sorted and sort it if it isn't
        echo(quiet, "\tVerifying if BED file is sorted")
        issorted = isbedsorted(bedutils, bed, noheader=noheader)

        if issorted is None:
            sys.stderr.write('Error in bedutils (issorted)\n')
            return 5
        elif issorted is False:
            echo(quiet, '\tSorting BED file')
            bed = os.path.join(tmpdir, os.path.basename(files[-1])[:-4] + '.sorted.bed')
            bed = sortbed(sort, files[-1], bed, ncores=ncores, maxmem=maxmem, tmpdir=tmpdir)

            if bed is None:
                sys.stderr.write('Error while sorting BED file\n')
                return 6
            files.append(bed)

    if nodup:
        echo(quiet, '\tRemoving PCR duplicate reads')
        bed = os.path.join(tmpdir, os.path.basename(bed)[:-4] + '.nodup.bed')
    else:
        bed = None

    nureads = filter_unique_reads(bedutils, files[-1], dest=bed, noheader=noheader)
    if nureads is None:
        sys.stderr.write('Error in bedutils (unique)\n')
        return 7
    elif nodup:
        files.append(bed)

    return files, nreads, nureads, chrms, len_mean


def getbackground(nreads, genomesize, binsize, ci):
    """
    Compute the background threshold using a simple Poisson distribution
    :param nreads:      Number of reads
    :param genomesize:  Size of the genome
    :param binsize:     Size of each bin
    :param ci:          Confidence interval
    :return:            int
    """
    # Lambda = number of reads / (genome size * bins size)
    lambda_l = float(nreads) / genomesize * binsize

    for bg in range(int(lambda_l)+1, 10001):
        # We only care for value greater than lambda (otherwise, we would get the threshold before the Poisson peak)
        p = poisson.pmf(bg, lambda_l)

        if p < (1 - ci):
            return bg

    return None


def get_genomeinfo(genome, chrms=None):
    """
    Return a dictionary containing the size of each chromosome
    If a list of chromosomes (found in the BED file),
    the functions returns a dictionary of the considered (known) chromosomes and the path to a temporary genome file
    :param genome:  Genome assembly
    :param chrms:   List of chromosomes (optional)
    :return:        Tuple (dictionary, file path)
    """
    path = os.path.join(os.path.dirname(__file__), 'genomes', genome)
    genomeinfo = {}
    with open(path) as f:
        for line in f:
            chrm, size = line.rstrip().split('\t')
            genomeinfo[chrm] = int(size)

    if chrms is None:
        return genomeinfo, None

    if len([i for i in chrms if i in genomeinfo]):
        pass
    else:
        genomeinfo2 = {}  # Adapted chromosomes for current BED file
        for chrm in chrms:

            if chrm.isdigit():
                string = 'chr' + chrm
                if string in genomeinfo:
                    # "1" -> "chr1"
                    genomeinfo2[chrm] = genomeinfo[string]
                    continue

                string = 'Chr' + chrm
                if string in genomeinfo:
                    # "1" -> "Chr1"
                    genomeinfo2[chrm] = genomeinfo[string]
                    continue
            elif chrm.lower() in genomeinfo:
                # "Chr1" -> "chr1"
                genomeinfo2[chrm] = genomeinfo[chrm.lower()]
                continue
            elif chrm[-3:].lower() == '.fa' and chrm[:-3] in genomeinfo:
                # "chr1.fa" -> "chr1"
                genomeinfo2[chrm] = genomeinfo[chrm[:-3]]
                continue

            if chrm[:3].lower() == 'chr':
                string = 'Chr' + chrm[3:]
                if string in genomeinfo:
                    # "chr1" -> "Chr1"
                    genomeinfo2[chrm] = genomeinfo[string]
                    continue

            # Ref chroms that contains the first characters of <chrm>
            c = [ref for ref in genomeinfo.keys() if chrm.startswith(ref)]
            if not c:
                continue

            # The list can contains several elements: we select the longest one
            # (chr11||hg19 starts with "chr11" and "chr1")
            genomeinfo2[chrm] = genomeinfo[max(c)]

        if genomeinfo2:
            # We managed to convert at least one chromosome
            genomeinfo = genomeinfo2
        else:
            # Well, we tried
            genomeinfo = None

    if genomeinfo:
        fd, tmppath = tempfile.mkstemp()
        os.close(fd)

        with open(tmppath, 'w') as fh:
            for chrm, size in sorted(genomeinfo.items()):
                fh.write('{}\t{}\n'.format(chrm, size))

        return genomeinfo, tmppath
    else:
        return genomeinfo, None


def sortbed(sort, bed, dest, ncores=1, maxmem=0, tmpdir=None):
    """
    Sort a BED file by chromosome, start position, and strand
    :param sort:        Sort binary
    :param bed:         BED file
    :param dest:        Output file
    :param ncores:      Number of cores
    :param maxmem:      Maximum buffer
    :param tmpdir:      Temporary directory
    :return:            None if an error occurs, the output file otherwise
    """
    cmd = [sort, '-t\t', '-k1,1', '-k2,2n', '-k6,6', bed]

    if ncores > 1:
        cmd += ['--parallel={}'.format(ncores)]

    if maxmem > 0:
        cmd += ['--buffer-size={}G'.format(maxmem)]

    if tmpdir is not None:
        cmd += ['-T', tmpdir]

    with open(dest, 'w') as fh, open(os.devnull, 'w') as dn:
        pop = Popen(cmd, stdout=fh, stderr=dn)
        pop.wait()

    if os.path.isfile(dest):
        if pop.returncode == 0:
            return dest
        os.unlink(dest)
    return None
