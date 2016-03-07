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
    parser.add_argument('-f')
    parser.add_argument('-b')

    # model
    parser.add_argument('--model_type')
    parser.add_argument('--alpha')
    parser.add_argument('--bound_low')
    parser.add_argument('--bound_high')

    args = parser.parse_args()

    # create the working dir

    os.mkdir('inputs')
    os.mkdir('job_outputs')
    os.mkdir('galaxy_outputs')

    # STACKS_archive
    # check if zipped files are into the tab

    extract_compress_files(args.f, "inputs")

    # create the populations command input line

    cmd_line = ['rxstacks']
    cmd_line.extend(["-b", args.b])
    cmd_line.extend(["-P", "inputs"])
    #cmd_line.extend(["-o", "job_outputs"])
    cmd_line.extend(["-o", "galaxy_outputs"])
    cmd_line.extend(["--model_type", args.model_type])

    if args.model_type == "snp":
        cmd_line.extend(["--alpha", args.alpha])
    elif args.model_type == "bounded":
        cmd_line.extend(["--alpha", args.alpha])
        cmd_line.extend(["--bound_low", args.bound_low])
        cmd_line.extend(["--bound_high", args.bound_high])
    else:
        pass

    print "[CMD_LINE] "+' '.join(cmd_line)


    # command with dependencies installed
    p = subprocess.Popen(cmd_line, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)

    (stdoutput, stderror) = p.communicate()

    print stdoutput
    print stderror


if __name__ == '__main__':
    __main__()


			
