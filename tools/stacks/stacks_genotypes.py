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
    parser.add_argument('-P')
    parser.add_argument('-b')
    parser.add_argument('-c')
    parser.add_argument('-t')
    parser.add_argument('-o')
    parser.add_argument('-e')
    parser.add_argument('--active_advanced')
    parser.add_argument('-r')
    parser.add_argument('-m')
    parser.add_argument('-B')
    parser.add_argument('-W')
    parser.add_argument('--active_autocorrect')
    parser.add_argument('--min_hom_seqs')
    parser.add_argument('--min_het_seqs')
    parser.add_argument('--max_het_seqs')

    # multifile management

    parser.add_argument('--logfile')
    parser.add_argument('--compress_output')

    # additionnal outputs

    parser.add_argument('--total_output')

    options = parser.parse_args()

        # create the working dir

    os.mkdir('job_outputs')
    os.mkdir('galaxy_outputs')

    os.chdir('job_outputs')

    # edit the command line

    cmd_line = []
    cmd_line.append("genotypes")

    # STACKS_archive
    # check if zipped files are into the tab

    extract_compress_files(options.P, os.getcwd())

    # create the genotypes command input line

    cmd_line.extend(["-b", options.b, "-P", os.getcwd()])   

    # create the genotypes command line

    if options.e:
       cmd_line.extend(["-e", options.e])
    if options.c == 'true':
       cmd_line.append("-c")
    if options.t:
        cmd_line.extend(["-t", options.t])
    if options.o:
        cmd_line.extend(["-o", options.o])

    # if advanced is activate
    if options.active_advanced == "true":
        cmd_line.extend(["-r", options.r])
        cmd_line.extend(["-m", options.m])
        if options.B:
            cmd_line.extend(["-B", options.B])
        if options.W:
            cmd_line.extend(["-W", options.W])

    # if autocorrect is activate
    if options.active_autocorrect == "true":
        cmd_line.extend(["--min_hom_seqs", options.min_hom_seqs])
        cmd_line.extend(["--min_het_seqs", options.min_het_seqs])
        cmd_line.extend(["--max_het_seqs", options.max_het_seqs])

    # command with dependencies installed
    print "[CMD]:"+' '.join(cmd_line)
    subprocess.call(cmd_line)

    # postprocesses
    try:
        shutil.copy('batch_1.haplotypes_1.tsv', options.logfile)
    except:
        sys.stderr.write('Error in genotypes execution; Please read the additional output (stdout)\n')
        sys.exit(1)

    # copy all files inside tmp_dir into workdir

    list_files = glob.glob('*')

    # if compress output is total

    if options.compress_output == 'total':
        mytotalzipfile = zipfile.ZipFile('total.zip.temp', 'w')

        for i in list_files:
            if re.search('^batch', os.path.basename(i)) \
                and not re.search("\.tsv$", os.path.basename(i)) \
                or re.search(".*_[0-9]*\.tsv$", os.path.basename(i)) \
                or re.search('.*genotypes.*', os.path.basename(i)):
                mytotalzipfile.write(i, os.path.basename(i))

        # return the unique archive

        shutil.move('total.zip.temp', options.total_output)

    # if compress output is default
    if options.compress_output == 'default':
        for i in list_files:
            if re.search('^batch', os.path.basename(i)) \
                and not re.search("\.tsv$", os.path.basename(i)) \
                or re.search(".*_[0-9]*\.tsv$", os.path.basename(i)) \
                or re.search('.*genotypes.*', os.path.basename(i)):
                shutil.move(i, '../galaxy_outputs')


if __name__ == '__main__':
    __main__()

			
