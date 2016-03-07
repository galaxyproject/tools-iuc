#!/usr/bin/env python

import sys, re
import os
import tempfile
import shutil, subprocess, glob
import optparse
from os.path import basename
import zipfile, tarfile, gzip
from galaxy.datatypes.checkers import *
from stacks import *

"""

Created by Yvan Le Bras
yvan.le_bras@irisa.fr

Last modifications : 02/17/2014

WARNING :

STACKS_sort_read_pairs.py needs:

- STACKS scripts in your $PATH

These scripts are available after compiling the sources of STACKS :

http://creskolab.uoregon.edu/stacks/

or with the galaxy_stacks package in the Genouest toolshed


"""

def __main__():
	

	# create the working dir
	os.mkdir("sort_read_outputs")
	os.mkdir("assembly_outputs")
	os.mkdir("samples_inputs")
        os.mkdir("stacks_inputs")

	# arguments recuperation
	parser = optparse.OptionParser()
	parser.add_option("-a")
	parser.add_option("-e")
	parser.add_option("-b")
	parser.add_option("-c")
   	parser.add_option("-d")
   	parser.add_option("-o")
	(options, args) = parser.parse_args()

	# edit the command line
	cmd_line1 = ["sort_read_pairs.pl"]

	#parse config files and create symlink if individual files are selected

	# STACKS_archive
	# check if zipped files are into the tab
	extract_compress_files(options.a, os.getcwd()+"/stacks_inputs")

	# samples_archive
	# check if zipped files are into the tab and change tab content
	extract_compress_files(options.e, os.getcwd()+"/samples_inputs")
		
	# create the sort_read_pairs command input line
	cmd_line1.extend(["-p", "stacks_inputs", "-s", "samples_inputs", "-o", "sort_read_outputs"])

	if options.b:
		cmd_line1.extend(["-w", options.b])
	if options.c:
		cmd_line1.extend(["-r", options.c])

	# exec command line 1 
	p1 = subprocess.Popen(cmd_line1)
        p1.communicate()

	# parse all files list and remove empty files from the output dir
        for fasta_file in glob.glob("sort_read_outputs/*"):
                if os.stat(fasta_file).st_size==0:
                        print "File "+fasta_file+" is empty"
                        os.remove(fasta_file)
        

	# create the exec_velvet command input line
	cmd_line2 = ["exec_velvet.pl"]
	cmd_line2.extend(["-s", "sort_read_outputs", "-o", "assembly_outputs"])
        cmd_line2.append("-c")

	if options.d:
		cmd_line2.extend(["-M", options.d])

	# version
	#cmd = 'sort_read_pairs.pl'+cmd_files+" "+cmd_options+" 2>&1"
	#cmd2 = 'exec_velvet.pl'+cmd_files2+" -c -e /softs/local/velvet/velvet_1.2.03/ "+cmd_options2+" 2>&1"
	
	# launch the command line 2 
        p2 = subprocess.Popen(cmd_line2)
        p2.communicate()

	# get collated.fa file
        try:
	    shutil.copy("assembly_outputs/collated.fa", options.o)
	except:
            print "No result file"
            sys.exit(1)

if __name__ == "__main__": __main__()
