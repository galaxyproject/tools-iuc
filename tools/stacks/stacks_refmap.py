#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import re
import os
import tempfile
import shutil
import subprocess
import glob
import optparse
from os.path import basename
import zipfile
import tarfile
import gzip
from galaxy.datatypes.checkers import *
from stacks import *


def __main__():

    # arguments recuperation

    parser = optparse.OptionParser()
    parser.add_option('-p')
    parser.add_option('-r')
    parser.add_option('-s')
    parser.add_option('-O')
    parser.add_option('-n')
    parser.add_option('-m')
    parser.add_option('-T')
    parser.add_option('--bound_low')
    parser.add_option('--bound_high')
    parser.add_option('--alpha')
    parser.add_option('--logfile')
    parser.add_option('--compress_output')
    parser.add_option('--catalogsnps')
    parser.add_option('--catalogalleles')
    parser.add_option('--catalogtags')

    # additionnal outputs

    parser.add_option('--total_output')
    parser.add_option('--tags_output')
    parser.add_option('--snps_output')
    parser.add_option('--alleles_output')
    parser.add_option('--matches_output')
    (options, args) = parser.parse_args()

        # create working directories

    os.mkdir('inputs')
    os.mkdir('job_outputs')
    os.mkdir('galaxy_outputs')

    cmd_line = []
    cmd_line.append('ref_map.pl')

    # if genetic map

    if options.p:

        # parse config files

        tab_parent_files = galaxy_config_to_tabfiles_for_STACKS(options.p)

        # check if zipped files are into the tab and change tab content

        extract_compress_files_from_tabfiles(tab_parent_files, 'inputs')

        # check files extension (important to have .sam files)

        check_sam_extension_and_add(tab_parent_files, 'inputs')

        # create symlink into the temp dir

        create_symlinks_from_tabfiles(tab_parent_files, 'inputs')

        # create the command input line

        for key in tab_parent_files:

            # if is a file (skip repository created after a decompression)

            if os.path.isfile('inputs/'+key):
                cmd_line.extend(['-p', os.path.normpath('inputs/'+key)])

    # if genetic map with progeny files

    if options.r:

        # parse config files

        tab_progeny_files = galaxy_config_to_tabfiles_for_STACKS(options.r)

        # check if zipped files are into the tab and change tab content

        extract_compress_files_from_tabfiles(tab_progeny_files, 'inputs')

        # check files extension (important to have .sam files)

        check_sam_extension_and_add(tab_progeny_files, 'inputs')

        # create symlink into the temp dir

        create_symlinks_from_tabfiles(tab_progeny_files, 'inputs')

        for key in tab_progeny_files:

            # if is a file (skip repository created after a decompression)

            if os.path.isfile('inputs/' + key):
                cmd_line.extend(['-r', 'inputs/' + key])

    # parse config files and create symlink if individual files are selected

    if options.s:

        # parse config files

        tab_individual_files = galaxy_config_to_tabfiles_for_STACKS(options.s)

        # check if zipped files are into the tab and change tab content

        extract_compress_files_from_tabfiles(tab_individual_files, 'inputs')

        # check files extension (important to have .sam files)

        check_sam_extension_and_add(tab_individual_files, 'inputs')

        # create symlink into the temp dir

        create_symlinks_from_tabfiles(tab_individual_files, 'inputs')

        # create the command input line

        for key in tab_individual_files:
            cmd_line.extend(['-s', 'inputs/' + key])

    # create the options command line

    cmd_line.extend([
        '-S',
        '-b', '1',
        '-T', '4',
        '-o', 'job_outputs',
        '-n', options.n,
        '-m', options.m,
        '-T', options.T,
        ])

    if options.O:
        cmd_line.extend(['-O', options.O])

    if options.bound_low:
        cmd_line.extend(['--bound_low', options.bound_low])

    if options.bound_high:
        cmd_line.extend(['--bound_high', options.bound_high])

    if options.alpha:
        cmd_line.extend(['--alpha', options.alpha])

    # execute job

    print '[COMMAND LINE]' + ' '.join(cmd_line)

    p = subprocess.Popen(cmd_line, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)

    (stdoutput, stderror) = p.communicate()

    print stdoutput
    print stderror

        # postprocesses

    try:
        shutil.move('job_outputs/ref_map.log', options.logfile)
    except:
        sys.stderr.write('Error in ref_map execution; Please read the additional output (stdout)\n')

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

    # copy all files inside tmp_dir into workdir

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

			
