#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import re
import tempfile
import subprocess
import glob
import shutil
import argparse
from os.path import basename
import zipfile
import tarfile
import gzip
from stacks import *


def __main__():

    # arguments recuperation

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_type')
    parser.add_argument('--input_enzyme')
    parser.add_argument('--input_single')
    parser.add_argument('--input_paired1')
    parser.add_argument('--input_paired2')
    parser.add_argument('--inputype')
    parser.add_argument('--sample_name')
    parser.add_argument('--barcode')
    parser.add_argument('--output_choice')
    parser.add_argument('--output_archive')
    parser.add_argument('--enzyme1')
    parser.add_argument('--enzyme2')
    parser.add_argument('--outype')
    parser.add_argument('--qualitenc')
    parser.add_argument('-D', action='store_true')
    parser.add_argument('--discard_file')
    parser.add_argument('-t')
    parser.add_argument('-q', action='store_true')
    parser.add_argument('--activate_advanced_options')
    parser.add_argument('-r', action='store_true')
    parser.add_argument('-w', default='0.15')
    parser.add_argument('-s', default='10')
    parser.add_argument('-c', action='store_true')
    parser.add_argument('--inline_null', action='store_true')
    parser.add_argument('--index_null', action='store_true')
    parser.add_argument('--inline_inline', action='store_true')
    parser.add_argument('--index_index', action='store_true')
    parser.add_argument('--inline_index', action='store_true')
    parser.add_argument('--index_inline', action='store_true')
    parser.add_argument('--logfile')
    options = parser.parse_args()

    # create the working dir
    #os.mkdir('inputs')
    os.mkdir('job_outputs')
    os.mkdir('galaxy_outputs')

    cmd_line = []
    cmd_line.append('process_radtags')
    #cmd_line.extend(['-p', 'inputs'])
    cmd_line.extend(['-i', options.inputype])
    cmd_line.extend(['-b', options.barcode])

    # parse config files and create symlink into the temp dir

    if options.input_type == 'single':

        # load the config file
        input_single = options.input_single

        # parse the input_file to extract filenames and filepaths
        tab_files = galaxy_config_to_tabfiles(input_single)

        # create symlink into the temp dir
        #create_symlinks_from_tabfiles(tab_files, 'inputs')
        cmd_line.extend(['-f', tab_files_paired.values()[0]])
    else:

        # load config files
        input_paired1 = options.input_paired1
        input_paired2 = options.input_paired2

        # parse the input_file to extract filenames and filepaths

        tab_files_paired1 = galaxy_config_to_tabfiles(input_paired1)
        tab_files_paired2 = galaxy_config_to_tabfiles(input_paired2)

        # create symlinks into the temp dir

        #create_symlinks_from_tabfiles(tab_files_paired1, 'inputs')
        #create_symlinks_from_tabfiles(tab_files_paired2, 'inputs')
        cmd_line.extend(['-1', tab_files_paired1.values()[0], '-2', tab_files_paired2.values()[0]])

        cmd_line.append('-P')

    # test nb enzyme
    if options.input_enzyme == '1':
        cmd_line.extend(['-e', options.enzyme1])

    if options.input_enzyme == '2':
        cmd_line.extend(['--renz_1', options.enzyme1, '--renz_2', options.enzyme2])

    # quality
    cmd_line.extend(['-E', options.qualitenc])

    # outputs
    cmd_line.extend(['-o', 'job_outputs/'])
    cmd_line.extend(['-y', options.outype])
    
    # test capture discards
    if options.D:
        cmd_line.append('-D')
    
    # optional options
    if options.activate_advanced_options == "true":

        if options.q:
            cmd_line.append('-q')
        if options.r:
            cmd_line.append('-r')
        
        cmd_line.extend(['-w', options.w, '-s', options.s])
        
        if options.c:
            cmd_line.append('-c')
        if options.t != '-1':
            cmd_line.extend(['-t', options.t])
        if options.inline_null:
            cmd_line.append('--inline_null')
        if options.index_null:
            cmd_line.append('--index_null')
        if options.inline_inline:
            cmd_line.append('--inline_inline')
        if options.index_index:
            cmd_line.append('--index_index')
        if options.inline_index:
            cmd_line.append('--inline_index')
        if options.index_inline:
            cmd_line.append('--index_inline')

    print '[CMD_LINE] : ' + str(cmd_line)
    print '[CMD_LINE] : ' + ' '.join(cmd_line)

    p = subprocess.call(cmd_line)

    # postprocesses

    try:
        shutil.move('job_outputs/process_radtags.log', options.logfile)
    except:
        sys.stderr.write('Error in process_radtags execution; Please read the additional output (stdout)\n')
        sys.exit(1)        

    if options.discard_file:
        discards_file_name = glob.glob('job_outputs/*.discards')[0]
        shutil.move(discards_file_name, options.discard_file)

    # manage outputs names

    change_outputs_procrad_name(os.getcwd() + '/job_outputs', options.sample_name)

    # generate additional output archive file

    if options.output_choice != '1':
        generate_additional_archive_file(os.getcwd() + '/job_outputs', options.output_archive)

    # if user has not choose the only zip archive

    if options.output_choice != '3':
        list_files = glob.glob('job_outputs/*')
        for i in list_files:
            shutil.move(i, 'galaxy_outputs')


if __name__ == '__main__':
    __main__()


