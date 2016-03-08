"""

STACKS METHODS FOR GALAXY

Created by Cyril Monjeaud & Yvan Le Bras
Cyril.Monjeaud@irisa.fr
yvan.le_bras@irisa.fr

Last modifications : 01/22/2014


"""

import os, sys, re, shutil
import glob 
import collections
import gzip, zipfile, tarfile
import subprocess
from galaxy.datatypes.checkers import *


"""

STACKS COMMON METHODS

galaxy_config_to_tabfiles(input_config)
galaxy_config_to_tabfiles_for_STACKS(input_config)
extract_compress_files_from_tabfiles(tab_files, tmp_input_dir)
create_symlinks_from_tabfiles(tab_files, tmp_input_dir)

"""
def galaxy_config_to_tabfiles(input_config):

	tab_files=collections.OrderedDict()
	for line in open(input_config, "r").readlines():
		if line.strip() != '':
			extract=line.strip().split("::")
			tab_files[extract[0].replace("(", ".").replace(" ", ".").replace(")", "").replace(":", ".").replace("/", ".")]=extract[1]

	# tabfiles[name]-> path
	return tab_files


def galaxy_config_to_tabfiles_for_STACKS(input_config):

	tab_files=collections.OrderedDict()
	for line in open(input_config, "r").readlines():
		if line.strip() != '':
			extract=line.strip().split("::")
			parse_name=re.search("STACKS.*\((.*\.[ATCG]*).*\)$", extract[0])
			# rename galaxy name in a short name
			if parse_name:				
				extract[0]=parse_name.groups(1)[0]

			tab_files[extract[0].replace("(", ".").replace(" ", ".").replace(")", "").replace(":", ".").replace("/", ".")]=extract[1]
			
	# tabfiles[name]-> path
	return tab_files


def extract_compress_files_from_tabfiles(tab_files, tmp_input_dir):
	
	# for each file
	for key in tab_files.keys():
		#test if is zip file
		if (check_zip( tab_files[key] )):

			# extract all files names and added it in the tab
			myarchive = zipfile.ZipFile(tab_files[key], 'r')
			for i in myarchive.namelist():
				tab_files[i]=tmp_input_dir+"/"+i

			# extract all files
			myarchive.extractall(tmp_input_dir)

			#remove compress file from the tab
			del tab_files[key]

		#test if is tar.gz file
		else:
			if tarfile.is_tarfile( tab_files[key] ) and check_gzip( tab_files[key] ):
				# extract all files names and added it in the tab
				mygzfile = tarfile.open(tab_files[key], 'r')
			
				for i in mygzfile.getnames():
					tab_files[i]=tmp_input_dir+"/"+i
		
				# extract all files
				mygzfile.extractall(tmp_input_dir)

				#remove compress file from the tab
				del tab_files[key]



def create_symlinks_from_tabfiles(tab_files, tmp_input_dir):

	for key in tab_files.keys():
		#print "file single: "+key+" -> "+tab_files[key]
		#create a sym_link in our temp dir
		if not os.path.exists(tmp_input_dir+'/'+key):
			os.symlink(tab_files[key], tmp_input_dir+'/'+key)


"""

PROCESS RADTAGS METHODS

generate_additional_file(tmp_output_dir, output_archive)

"""

def change_outputs_procrad_name(tmp_output_dir, sample_name):

	list_files = glob.glob(tmp_output_dir+'/*')

	for file in list_files:
		# change sample name
		new_file_name=os.path.basename(file.replace("_",".").replace("sample", sample_name))

		# transform .fa -> .fasta or .fq->.fastq
		if os.path.splitext(new_file_name)[1] == ".fa":
			new_file_name = os.path.splitext(new_file_name)[0]+'.fasta'
		else:
			new_file_name = os.path.splitext(new_file_name)[0]+'.fastq'

		shutil.move(tmp_output_dir+'/'+os.path.basename(file), tmp_output_dir+'/'+new_file_name)


def generate_additional_archive_file(tmp_output_dir, output_archive):

	list_files = glob.glob(tmp_output_dir+'/*')

        myzip=zipfile.ZipFile("archive.zip.temp", 'w', allowZip64=True)

	# for each fastq file
	for fastq_file in list_files:
		# add file to the archive output
		myzip.write(fastq_file, os.path.basename(fastq_file))

	shutil.move("archive.zip.temp", output_archive)


"""

DENOVOMAP METHODS

check_fastq_extension_and_add(tab_files, tmp_input_dir)

"""

def check_fastq_extension_and_add(tab_files, tmp_input_dir):
	
	# for each file
	for key in tab_files.keys():

		if not re.search("\.fq$", key) and not re.search("\.fastq$", key) and not re.search("\.fa$", key) and not re.search("\.fasta$", key):
			# open the file
			myfastxfile=open(tab_files[key], 'r')
			
			# get the header
			line = myfastxfile.readline()
		        line = line.strip()

			# fasta rapid test
			if line.startswith( '>' ):
                        	tab_files[key+".fasta"]=tab_files[key]
				del tab_files[key]
			# fastq rapid test
			elif line.startswith( '@' ):
				tab_files[key+".fq"]=tab_files[key]
				del tab_files[key]
			else:
				print "[WARNING] : your input file "+key+" was not extension and is not recognize as a Fasta or Fastq file"
		
			myfastxfile.close()


"""

REFMAP METHODS

"""

def check_sam_extension_and_add(tab_files, tmp_input_dir):
	
	# for each file
	for key in tab_files.keys():

		if not re.search("\.sam$", key):
			# add the extension
			tab_files[key+".sam"]=tab_files[key]
			del tab_files[key]






"""

PREPARE POPULATION MAP METHODS

generate_popmap_for_denovo(tab_files, infos_file, pop_map)
generate_popmap_for_refmap(tab_fq_files, tab_sam_files, infos_file, pop_map)


"""
def generate_popmap_for_denovo(tab_files, infos_file, pop_map):

	# initiate the dict : barcode -> tab[seq]
	fq_name_for_barcode={}

	for key in tab_files:
		single_barcode=re.search("([ATCG]*)(\.fq|\.fastq)", key).groups(0)[0]
		fq_name_for_barcode[single_barcode]=key

	# open the infos file and output file
	my_open_info_file=open(infos_file, 'r')
	my_output_file=open(pop_map, 'w')

	# conversion tab for population to integer
	pop_to_int=[]

	# write infos into the final output
	for line in my_open_info_file:

		parse_line=re.search("(^[ATCG]+)\t(.*)", line.strip())

		if not parse_line:
			print "[WARNING] Wrong input infos file structure : "+line
		else:
			barcode=parse_line.groups(1)[0]
			population_name=parse_line.groups(1)[1]

			# if its the first meet with the population
			if population_name not in pop_to_int:
				pop_to_int.append(population_name)		

			# manage ext if present, because the population map file should not have the ext
			if re.search("(\.fq$|\.fastq$)", fq_name_for_barcode[barcode]):
				fqfile=os.path.splitext(fq_name_for_barcode[barcode])[0]
			else:
				fqfile=fq_name_for_barcode[barcode]


			# write in the file
			my_output_file.write(fqfile+"\t"+str(pop_to_int.index(population_name))+"\n")

	# close files
	my_output_file.close()
	my_open_info_file.close()




def generate_popmap_for_refmap(tab_fq_files, tab_sam_files, infos_file, pop_map):

	# initiate the dict : barcode -> tab[seq]
	seq_id_for_barcode={}

	# initiate the dict : barcode -> sam_name
	sam_name_for_barcode={}

	### Parse fastqfiles ###
	# insert my barcode into a tab with sequences ID associated
	for fastq_file in tab_fq_files.keys():
		single_barcode=re.search("([ATCG]*)(\.fq|\.fastq)", fastq_file).groups(0)[0]
        
        	# open the fasq file
        	open_fastqfile=open(tab_fq_files[fastq_file], 'r')
		
		# for each line, get the seq ID
		tab_seq_id=[]
		for line in open_fastqfile:
		    my_match_seqID=re.search("^@([A-Z0-9]+\.[0-9]+)\s.*", line)
		    if my_match_seqID:    
		        tab_seq_id.append(my_match_seqID.groups(0)[0])

		# push in a dict the tab of seqID for the current barcode
		seq_id_for_barcode[single_barcode]=tab_seq_id


	### Parse samfiles and get the first seq id ###
	for sam_file in tab_sam_files.keys():

		# open the sam file
        	open_samfile=open(tab_sam_files[sam_file], 'r')
		
		# get the first seq id
		first_seq_id=''
		for line in open_samfile:
			if not re.search("^@", line):
				first_seq_id=line.split("\t")[0]
				break


		# map with seq_id_for_barcode structure
		for barcode in seq_id_for_barcode:
			for seq in seq_id_for_barcode[barcode]:
				if seq == first_seq_id:
					#print "sam -> "+sam_file+" seq -> "+first_seq_id+" barcode -> "+barcode
					sam_name_for_barcode[barcode]=sam_file
					break

	# open the infos file and output file
	my_open_info_file=open(infos_file, 'r')
	my_output_file=open(pop_map, 'w')

	# conversion tab for population to integer
	pop_to_int=[]

	# write infos into the final output
	for line in my_open_info_file:
		parse_line=re.search("(^[ATCG]+)\t(.*)", line)

		if not parse_line:
			print "[WARNING] Wrong input infos file structure : "+line
		else:

			# if its the first meet with the population
			if parse_line.groups(1)[1] not in pop_to_int:
				pop_to_int.append(parse_line.groups(1)[1])		

			# manage ext if present, because the population map file should not have the ext
			if re.search("\.sam", sam_name_for_barcode[parse_line.groups(1)[0]]):
				samfile=os.path.splitext(sam_name_for_barcode[parse_line.groups(1)[0]])[0]
			else:
				samfile=sam_name_for_barcode[parse_line.groups(1)[0]]

			# write in the file
			my_output_file.write(samfile+"\t"+str(pop_to_int.index(parse_line.groups(1)[1]))+"\n")

	# close files
	my_output_file.close()
	my_open_info_file.close()


"""

STACKS POPULATION


"""


def extract_compress_files(myfile, tmp_input_dir):

        #test if is zip file
        if (check_zip( myfile )):

                # extract all files names and added it in the tab
                myarchive = zipfile.ZipFile(myfile, 'r')

                # extract all files
                myarchive.extractall(tmp_input_dir)


        #test if is tar.gz file
        else:
                # extract all files names and added it in the tab
                mygzfile = tarfile.open(myfile, 'r')

                # extract all files
                mygzfile.extractall(tmp_input_dir)


