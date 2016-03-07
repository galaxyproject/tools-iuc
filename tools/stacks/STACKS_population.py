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
    parser.add_argument('-M')
    parser.add_argument('-b')
    parser.add_argument('--vcf', action='store_true')
    parser.add_argument('--genepop', action='store_true')
    parser.add_argument('--structure', action='store_true')
    parser.add_argument('-e')
    parser.add_argument('--genomic', action='store_true')
    parser.add_argument('--fasta', action='store_true')
    parser.add_argument('--phase', action='store_true')
    parser.add_argument('--beagle', action='store_true')
    parser.add_argument('--plink', action='store_true')
    parser.add_argument('--phylip', action='store_true')
    parser.add_argument('--phylip_var', action='store_true')
    parser.add_argument('--write_single_snp', action='store_true')
    parser.add_argument('-k', action='store_true')

    # advanced options
    parser.add_argument('--advanced_options_activate')
    parser.add_argument('-B')
    parser.add_argument('-W')
    parser.add_argument('-r')
    parser.add_argument('-p')
    parser.add_argument('-m')
    parser.add_argument('-a')
    parser.add_argument('-f')
    parser.add_argument('-t')
    parser.add_argument('--p_value_cutoff')
    parser.add_argument('--window_size')
    parser.add_argument('--bootstrap')
    parser.add_argument('--bootstrap_reps')

    # multifile management
    parser.add_argument('--logfile')

    # outputs
    parser.add_argument('--ss')
    parser.add_argument('--s')

    # optional outputs
    parser.add_argument('--ov')
    parser.add_argument('--op')
    parser.add_argument('--ol')
    parser.add_argument('--of')
    parser.add_argument('--os')
    parser.add_argument('--oe')
    parser.add_argument('--om')
    parser.add_argument('--og') 

    parser.add_argument('--unphased_output')
    parser.add_argument('--markers_output')
    parser.add_argument('--phase_output')
    parser.add_argument('--fst_output')

    options = parser.parse_args()

    # create the working dir
    os.mkdir('job_outputs')
    os.mkdir('galaxy_outputs')

    os.chdir('job_outputs')

    # STACKS_archive
    # check if zipped files are into the tab
    extract_compress_files(options.P, os.getcwd())

    # create the populations command input line
    cmd_line=['populations']
    cmd_line.extend(['-b', options.b, '-P', os.getcwd(), '-M', options.M])

    if options.e:
        cmd_line.extend(['-e', options.e, options.genomic])

    # output options
    if options.vcf:
        cmd_line.append('--vcf')
    if options.genepop:
        cmd_line.append('--genepop')
    if options.structure:
        cmd_line.append('--structure')
    if options.fasta:
        cmd_line.append('--fasta')
    if options.phase:
        cmd_line.append('--phase')
    if options.beagle:
        cmd_line.append('--beagle')
    if options.plink:
        cmd_line.append('--plink')
    if options.phylip:
        cmd_line.append('--phylip')
    if options.phylip_var and options.phylip:
        cmd_line.append('--phylip_var')
    if options.write_single_snp and (options.genepop or options.structure):
        cmd_line.append('--write_single_snp')

    if options.k:
        cmd_line.extend(['-k', '--window_size', options.window_size])
   
    if options.advanced_options_activate == 'true':
        if options.B:
            cmd_line.extend(['-B', options.B])
        if options.W:
            cmd_line.extend(['-W', options.W])

        cmd_line.extend(['-r', options.r])
        cmd_line.extend(['-p', options.p])
        cmd_line.extend(['-m', options.m])
        cmd_line.extend(['-a', options.a])

    if options.f:
        cmd_line.extend(['-f', options.f, '--p_value_cutoff', options.p_value_cutoff])
    if options.bootstrap:
        cmd_line.extend(['--bootstrap', options.bootstrap, '--bootstrap_reps', options.bootstrap_reps])

    if options.t:
        cmd_line.extend(['-t', options.t])

    print "[CMD]:"+' '.join(cmd_line)
    subprocess.call(cmd_line)

    # postprocesses
    try:
        shutil.copy('batch_1.populations.log', options.logfile)
    except:
        sys.stderr.write('Error in population execution; Please read the additional output (stdout)\n')
        sys.exit(1)

    try:
        shutil.move(glob.glob('*.sumstats_summary.tsv')[0], options.ss)
    except:
        print "No sumstats summary file"

    try:
        shutil.move(glob.glob('*.sumstats.tsv')[0], options.s)
    except:
        print "No sumstats file"

    # move additionnal output files
    if options.vcf:
        try:
            shutil.move(glob.glob('*.vcf')[0], options.ov)
        except:
            print "No VCF files"

    if options.phylip:
        try:
            shutil.move(glob.glob('*.phylip')[0], options.op)
            shutil.move(glob.glob('*.phylip.log')[0], options.ol)
        except:
            print "No phylip file"

    if options.fasta:
        try:
            shutil.move(glob.glob('*.fa')[0], options.of)
        except:
            print "No fasta files"

    if options.structure:
        try:
            shutil.move(glob.glob('*.structure.tsv')[0], options.os)
        except:
            print "No structure file"

    if options.plink :
        try:
            shutil.move(glob.glob('*.ped')[0], options.oe)
            shutil.move(glob.glob('*.map')[0], options.om)
        except:
            print "No ped and map file"

    if options.genepop :
        try:
            shutil.move(glob.glob('*.genepop')[0], options.og)
        except:
            print "No genepop file"

    # copy all files inside tmp_dir into workdir or into an archive....
    list_files = glob.glob('*')

    markerszip = zipfile.ZipFile('markers.zip.temp', 'w',
                                 allowZip64=True)
    phasezip = zipfile.ZipFile('phase.zip.temp', 'w', allowZip64=True)
    unphasedzip = zipfile.ZipFile('unphased.zip.temp', 'w',
                                  allowZip64=True)
    fstzip = zipfile.ZipFile('fst.zip.temp', 'w', allowZip64=True)

    for i in list_files:
        # for each type of files
        if re.search("\.markers$", i):
            markerszip.write(i)
        elif re.search("phase\.inp$", i):
            phasezip.write(i)
        elif re.search("unphased\.bgl$", i):
            unphasedzip.write(i)
        elif re.search('fst', i):
            fstzip.write(i)
        else:
        # else return original files
            if re.search('^batch', os.path.basename(i)) \
                and not re.search("\.tsv$", os.path.basename(i)) \
                or re.search(".*_[0-9]*\.tsv$", os.path.basename(i)):
                shutil.move(i, '../galaxy_outputs')

    # close zip files
    markerszip.close()
    phasezip.close()
    unphasedzip.close()
    fstzip.close()

    # return archives
    shutil.move('fst.zip.temp', options.fst_output)
    if options.beagle:
        shutil.move('markers.zip.temp', options.markers_output)
        shutil.move('unphased.zip.temp', options.unphased_output)
    if options.phase:
        shutil.move('phase.zip.temp', options.phase_output)


if __name__ == '__main__':
    __main__()


			
