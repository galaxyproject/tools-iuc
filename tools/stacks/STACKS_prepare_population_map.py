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

Created by Cyril Monjeaud
Cyril.Monjeaud@irisa.fr

Last modifications : 01/10/2014

WARNING :

STACKS_denovomap.py needs:

- STACKS scripts in your $PATH

These scripts are available after compiling the sources of STACKS :

http://creskolab.uoregon.edu/stacks/

or with the galaxy_stacks package in the Genouest toolshed (http://toolshed.genouest.org)

"""
def __main__():
	
	# arguments recuperation
	parser = optparse.OptionParser()
	parser.add_option("-f")
	parser.add_option("-s")
	parser.add_option("-t")
	parser.add_option("-o")
	parser.add_option("-d")
	(options, args) = parser.parse_args()

	# create the working dir
	tmp_dir = tempfile.mkdtemp(dir=options.d)

	print tmp_dir
	#os.chdir(tmp_dir)

	# parse config files
	tab_fq_files=galaxy_config_to_tabfiles_for_STACKS(options.f)

	# check if zipped files are into the tab and change tab content
	extract_compress_files_from_tabfiles(tab_fq_files, tmp_dir)

	# generate population map for denovo map
	if not options.s:
		generate_popmap_for_denovo(tab_fq_files, options.t, options.o)
	else:
		# parse config files
		tab_sam_files=galaxy_config_to_tabfiles_for_STACKS(options.s)
		extract_compress_files_from_tabfiles(tab_sam_files, tmp_dir)
		generate_popmap_for_refmap(tab_fq_files, tab_sam_files, options.t, options.o)
	

	#clean up temp files
	shutil.rmtree( tmp_dir )

	
	



if __name__ == "__main__": __main__()
