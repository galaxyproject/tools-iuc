#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import re
import os
import tempfile
import shutil
import subprocess
import glob
import argparse
from os.path import basename
import zipfile
import tarfile
import gzip
from galaxy.datatypes.checkers import *
from stacks import *


def __main__():

    # arguments recuperation

    parser = argparse.ArgumentParser()
    parser.add_argument('-p')
    parser.add_argument('-b')
    parser.add_argument('-r')
    parser.add_argument('-s')
    parser.add_argument('-O')
    parser.add_argument('-m')
    parser.add_argument('-P')
    parser.add_argument('-M')
    parser.add_argument('-N')
    parser.add_argument('-n')
    parser.add_argument('-T')
    parser.add_argument('-t')
    parser.add_argument('-H')
    parser.add_argument('--bound_low')
    parser.add_argument('--bound_high')
    parser.add_argument('--alpha')
    parser.add_argument('--logfile')
    parser.add_argument('--compress_output')
    parser.add_argument('--catalogsnps')
    parser.add_argument('--catalogalleles')
    parser.add_argument('--catalogtags')

    # additionnal outputs
    parser.add_argument('--total_output')
    parser.add_argument('--tags_output')
    parser.add_argument('--snps_output')
    parser.add_argument('--alleles_output')
    parser.add_argument('--matches_output')

    options = parser.parse_args()

    # create working directories

    os.mkdir('inputs')
    os.mkdir('job_outputs')
    os.mkdir('galaxy_outputs')

    cmd_line = []
    cmd_line.append('denovo_map.pl')

    # if genetic map

    if options.p:

        # parse config files

        tab_parent_files = galaxy_config_to_tabfiles_for_STACKS(options.p)

        # check if zipped files are into the tab and change tab content

        extract_compress_files_from_tabfiles(tab_parent_files, 'inputs')

        # check files extension (important to have .fq or .fasta files)

        check_fastq_extension_and_add(tab_parent_files, 'inputs')

        # create symlink into the temp dir

        create_symlinks_from_tabfiles(tab_parent_files, 'inputs')

        # parse the input dir and store all file names into a tab

        fastq_files = []
        for fastq_file in glob.glob('inputs/*'):
            # if is a file (skip repository created after a decompression)
            if os.path.isfile(fastq_file):
                fastq_files.append(fastq_file)

        fastq_files.sort()

        # test if fastq are paired-end
        if options.b == 'true':
            for n in range(0, len(fastq_files), 2):
                cmd_line.extend(['-p', fastq_files[n]])
        else:
            for myfastqfile in fastq_files:
                cmd_line.extend(['-p', myfastqfile])

    # if genetic map with progeny files

    if options.r:

        # parse config files
        tab_progeny_files = galaxy_config_to_tabfiles_for_STACKS(options.r)

        # check if zipped files are into the tab and change tab content
        extract_compress_files_from_tabfiles(tab_progeny_files, 'inputs')

        # check files extension (important to have .fq or .fasta files)
        check_fastq_extension_and_add(tab_progeny_files, 'inputs')

        # create symlink into the temp dir
        create_symlinks_from_tabfiles(tab_progeny_files, 'inputs')

        for key in tab_progeny_files:

            # if is a file (skip repository created after a decompression)

            if os.path.isfile('inputs/' + key):
                cmd_line.extend(['-r', 'inputs/' + key])

    # if population is checked
    if options.s:

        tab_individual_files = galaxy_config_to_tabfiles_for_STACKS(options.s)

        # check if zipped files are into the tab and change tab content
        extract_compress_files_from_tabfiles(tab_individual_files, 'inputs')

        # check files extension (important to have .fq or .fasta files)
        check_fastq_extension_and_add(tab_individual_files, 'inputs')

        # create symlink into the temp dir
        create_symlinks_from_tabfiles(tab_individual_files, 'inputs')

        # create the command input line
        for key in tab_individual_files:

            # if is a file (skip repository created after a decompression)
            if os.path.isfile('inputs/' + key):
                cmd_line.extend(['-s', 'inputs/' + key])

    # create the command line
    cmd_line.extend([
        '-S',
        '-b',
        '1',
        '-T',
        '4',
        '-o',
        'job_outputs/'
        ])

    if options.O:
        cmd_line.extend(['-O', options.O])

    if options.m and options.m != '-1':
        cmd_line.extend(['-m', options.m])

    if options.P and options.P != '-1':
        cmd_line.extend(['-P', options.P])

    if options.M and options.M != '-1':
        cmd_line.extend(['-M', options.M])

    if options.N and options.N != '-1':
        cmd_line.extend(['-N', options.N])

    if options.n and options.n != '-1':
        cmd_line.extend(['-n', options.n])

    if options.T and options.T != '-1':
        cmd_line.extend(['-T', options.T])

    if options.t and options.t == 'true':
        cmd_line.append('-t')

    if options.H and options.H == 'true':
        cmd_line.append('-H')

    ## SNP model 
    if options.bound_low:
        cmd_line.extend(['--bound_low', options.bound_low])
        cmd_line.extend(['--bound_high', options.bound_high])

    if options.alpha:
        cmd_line.extend(['--alpha', options.alpha])

    # launch the command line
    print "[CMD_LINE] : "+' '.join(cmd_line)    

    p = subprocess.call(cmd_line)

    # postprocesses
    try:
        shutil.move('job_outputs/denovo_map.log', options.logfile)
    except:
        sys.stderr.write('Error in denovo_map execution; Please read the additional output (stdout)\n')
        sys.exit(1)

    # go inside the outputs dir
    os.chdir('job_outputs')

    # move files
    for i in glob.glob('*'):
        if re.search('catalog.snps.tsv$', i):
            shutil.copy(i, options.catalogsnps)
        if re.search('catalog.alleles.tsv$', i):
            shutil.copy(i, options.catalogalleles)
        if re.search('catalog.tags.tsv$', i):
            shutil.copy(i, options.catalogtags)

    list_files = glob.glob('*')

    # if compress output is total
    if options.compress_output == 'total':

        mytotalzipfile = zipfile.ZipFile('total.zip.temp', 'w',
                allowZip64=True)

        for i in list_files:
            mytotalzipfile.write(os.path.basename(i))

        # return the unique archive
        shutil.move('total.zip.temp', options.total_output)
    elif options.compress_output == 'categories':

    # if compress output is by categories
        mytagszip = zipfile.ZipFile('tags.zip.temp', 'w', allowZip64=True)
        mysnpszip = zipfile.ZipFile('snps.zip.temp', 'w', allowZip64=True)
        myalleleszip = zipfile.ZipFile('alleles.zip.temp', 'w', allowZip64=True)
        mymatcheszip = zipfile.ZipFile('matches.zip.temp', 'w', allowZip64=True)

        for i in list_files:
            # for each type of files
            if re.search("tags\.tsv$", i) and not re.search('batch', i):
                mytagszip.write(os.path.basename(i))
                os.remove(i)
            elif re.search("snps\.tsv$", i) and not re.search('batch', i):
                mysnpszip.write(os.path.basename(i))
                os.remove(i)
            elif re.search("alleles\.tsv$", i) and not re.search('batch', i):
                myalleleszip.write(os.path.basename(i))
                os.remove(i)
            elif re.search("matches\.tsv$", i) and not re.search('batch', i):
                mymatcheszip.write(os.path.basename(i))
                os.remove(i)
            else:
                shutil.move(os.path.basename(i), '../galaxy_outputs')

        # return archives....
        shutil.move('tags.zip.temp', options.tags_output)
        shutil.move('snps.zip.temp', options.snps_output)
        shutil.move('alleles.zip.temp', options.alleles_output)
        shutil.move('matches.zip.temp', options.matches_output)
    else:
    # else no compression
        for i in list_files:
            shutil.move(os.path.basename(i), '../galaxy_outputs')


if __name__ == '__main__':
    __main__()

			
