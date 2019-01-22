#!/usr/bin/env python

#--------------------------------------
#
#        snippy_core_wrapper.py
#
# This is an intermediary script between snippy-core.xml and snippy-core
# It:
#   - Copys the supplied zipped snippy output files to the working dir
#   - Untars them to their datafile name
#   - Builds the snippy-core command
#   - Runs the snippy-core command
#
#--------------------------------------

import os
import sys
import argparse
import subprocess
from shutil import copyfile

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--ref', help='Reference fasta', required=True)
    parser.add_argument('-i', '--indirs', help='Comma-separated list of input datasets')
    args = parser.parse_args()

    snippy_core_command_line = ['snippy-core', '--ref', args.ref]

    for indir in args.indirs.split(','):
        base_name = os.path.basename(indir)
        sys.stderr.write(indir + '\n')
        copyfile(indir, base_name)
        sys.stderr.write('\n'.join([indir, base_name]) + '\n')
        subprocess.Popen(['tar', '-xf', base_name])

    extracted_dirs = [f.path for f in os.listdir() if f.is_dir() ]

    for extracted_dir in extracted_dirs:
        snippy_core_command_line.append(extracted_dir)

    print(snippy_core_command_line)
    subprocess.Popen(snippy_core_command_line)

if __name__ == '__main__':
    main()
