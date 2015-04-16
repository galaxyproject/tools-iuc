import pprint
import os.path
import itertools
import re
import subprocess
import tempfile


def xmfa_parse(xmfa):
	start_regex = re.compile('> (?P<idx>[1-9]+):(?P<start>[0-9]+)-(?P<end>[0-9]+)')
	seq_info_list = {}
	with open(xmfa, 'r') as handle:
		seqidx = None
		seq_start = None

		for line in handle.readlines():
			match = start_regex.match(line)
			if match is not None:
				seqidx = match.group('idx')
				if seqidx not in seq_info_list:
					seq_info_list[seqidx] = {}

				seq_start = int(start_regex.match(line).group('start'))
				continue

			# = separates blocks
			if line.strip() == "=":
				continue

			# This should ONLY be sequence. We've skipped header, we've passed
			# >s and =s
			#
			# Checking that seqidx is not None ensure that we're within a
			# sequence block
			if not line.startswith('#') and seqidx is not None:
				if seq_start not in seq_info_list[seqidx]:
					seq_info_list[seqidx][seq_start] = ''
				seq_info_list[seqidx][seq_start] += line.strip()

	return_list = []
	for seqidx in sorted(seq_info_list.keys()):
		return_list.append(''.join([seq_info_list[seqidx][x] for x in sorted(seq_info_list[seqidx].keys())]))
	return return_list


def links(backbone_file, link_output):
	header_parsed = False
	links = {}
	with open(backbone_file, 'r') as handle:
		for line in handle:
			if not header_parsed:
				# There will be 2N where N is the number of genomes in our header_parsed
				header_parsed = line.split('\t')
				genome_count = len(header_parsed) / 2
				# Data structure to hold the links
				# A from-to list (from col/to col) which will be used to check links between genomes
				# https://docs.python.org/2/library/itertools.html
				# from_to for a genome_count of 3 looks like: [(0,  1),  (0,  2),  (1,  2)]
				# Preferred as it will expand to ANY number of comparisons
				from_to = list(itertools.combinations(range(genome_count), 2))
				header_parsed = True
			else:
				# Now down to actual parsing
				link_data = line.split()
				# Iterate over pairs of columns
				for x, y in from_to:
					# Access the link data in that column
					a_left = link_data[2 * x].strip()
					a_right = link_data[2 * x + 1].strip()
					b_left = link_data[2 * y].strip()
					b_right = link_data[2 * y + 1].strip()
					# if any of them are zero,  then we can continue,  as this
					# isn't a "true" link. Circos will plot links from "0 0",  so
					# any links with "0 0" need to be removed.
					if a_left != 0 and b_left != 0:
						key = '%s-%s' % (x, y)
						if key not in links:
							links[key] = []

						links['%s-%s' % (x, y)].append(
							map(int, [a_left, a_right, b_left, b_right]))
	return links

def reverse_complement(sequence):
	DNA_pairing_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G', '-':'-'}
	new_seq = ''
	for letter in sequence[::-1]:
		new_seq+=DNA_pairing_dict[letter]
	return new_seq


def percent_sequence_identity(seq1, seq2):
	if len(seq1) == 0 or len(seq2) == 0:
		return 0

	results = [seq1[i] == seq2[i]
				for i in range(min(len(seq1), len(seq2)))]
	# Count the trues, return that as % of length
	# Currently will underestimate sequences
	return float(results.count(True)) / len(results)

def write_link_file(names_list,  links,  link_output, directory):
	with open(os.path.join(directory, link_output),  'w') as handle:
		for key in links:
			key_from,  key_to = key.split('-') # Rememver,  we keyed on 1-3 0-2 3-2 / etc
			for link in links[key]: # List of from/to
		# Create the list by re-arranging the link_data
		# Circos will plot links from "0 0",  so any links with "0 0" need to be removed.
				if len(link) == 5:
					if int(link[0]) > 0 and int(link[2])>0:
						data = [names_list[int(key_from)]] + link[0:2] + [names_list[int(key_to)]] + link[2:4]
						if float(link[4])<0.3:
							handle.write(' '.join(data) + " color=lblue \n")
							#print ' '.join(data) + " color=lblue \n"
						elif float(link[4])<0.6:
							handle.write(' '.join(data) + " color=blue \n")
							#print ' '.join(data) + " color=lblue \n"
						else:
							handle.write(' '.join(data) + " color=dblue \n")
							#print ' '.join(data) + " color=lblue \n"
					elif int(link[0]) < 0:
						link[0] = str(-1*int(link[0]))
						link[1] = str(-1*int(link[1]))
						data = [names_list[int(key_from)]] + link[0:2] + [names_list[int(key_to)]] + link[2:4]
						if float(link[4])<0.3:
							handle.write(' '.join(data) + " color=lred \n")
							#print ' '.join(data) + " color=lblue \n"
						elif float(link[4])<0.6:
							handle.write(' '.join(data) + " color=red \n")
							#print ' '.join(data) + " color=lblue \n"
						else:
							handle.write(' '.join(data) + " color=dred \n")
							#print ' '.join(data) + " color=lblue \n"
					else:
						link[2] = str(-1*int(link[2]))
						link[3] = str(-1*int(link[3]))
						data = [names_list[int(key_from)]] + link[0:2] + [names_list[int(key_to)]] + link[2:4]
						if float(link[4])<0.3:
							handle.write(' '.join(data) + " color=lred \n")
							#print ' '.join(data) + " color=lblue \n"
						elif float(link[4])<0.6:
							handle.write(' '.join(data) + " color=red \n")
							#print ' '.join(data) + " color=lblue \n"
						else:
							handle.write(' '.join(data) + " color=dred \n")
							#print ' '.join(data) + " color=lblue \n"



def karyotype(seq_file, karyotype_name, directory):
	i=0
	colors_list = ['red', 'blue', 'green', 'orange', 'violet', 'brown']
	genome_list = []
	names_list = []
	seq_file = open(seq_file, 'r')
	preceding_line_was_header = 0
	for line in seq_file.readlines():
		if preceding_line_was_header == 1:
			genome_tuple = (name, len(line))
			genome_list += [genome_tuple]
		preceding_line_was_header = 0
		name = ''
		genome_tuple = ('', '')
		if line[0] == '>':
			preceding_line_was_header = 1
		#Circos does not like ">" in the genome names.
		name = line[1:].strip()
	seq_file.close()
	karyotype_file = open(os.path.join(directory, karyotype_name), 'w')
	for tuple in genome_list:
		names_list+=[tuple[0]]
		karyotype_file.write('chr - '+(tuple[0]+' ')*2+str(0)+' '+str(tuple[1]-1)+ ' ' + str(colors_list[i]) + ' \n')
		i+=1
	karyotype_file.close()
	return names_list

def add_pct_identity(link_dict, sequence_file, alignment_list):
	with open(sequence_file, 'r') as handle:
		lines_list = handle.readlines()
		for key in link_dict:
			list_of_links = link_dict[key]
			key_pair = map(int, key.split('-'))

			for element in list_of_links:
				first_index_set = [
					min(abs(element[0]), abs(element[1])),
					max(abs(element[0]), abs(element[1]))]
				second_index_set = [
					min(abs(element[2]), abs(element[3])),
					max(abs(element[2]), abs(element[3]))]

				sect_1 = alignment_list[key_pair[0]][first_index_set[0]:first_index_set[1]]
				sect_2 = alignment_list[key_pair[1]][second_index_set[0]:second_index_set[1]]
				print key_pair[0], first_index_set[0], first_index_set[1], sect_1
				print key_pair[0], second_index_set[0], second_index_set[1], sect_2
				print

				if element[0] < 0:
					sect_1 = reverse_complement(sect_1)
				if element[2] < 0:
						sect_2 = reverse_complement(sect_2)

				if first_index_set[0] != 0 and second_index_set != 0:
					element += [percent_sequence_identity(sect_1, sect_2)]
					#pprint.pprint(element)
			#for item in key_pair:
			#	print lines_list[2*int(item)+1].strip()[]
		pprint.pprint(link_dict)
		return link_dict


def run_progressiveMauve(sequence_file):
	# Tempfiles will need to be cleaned up when they're finished with, but for testing we'll just leave it here.
	#tmp = tempfile.NamedTemporaryFile(delete=False)
	# Then use tmp.name to access the file path
	#
	tmp = "output.xmfa"
	subprocess.check_call(['progressiveMauve', '--output=%s' % tmp, sequence_file])


if __name__ == '__main__':
	directory = 'C:\\Users\\User\\Desktop\\491 Scripts\\test_2\\'
	destination_directory = 'C:\\Users\\User\\Desktop\\491 Scripts\\test_2\\'
	image_direct = 'C:\\Users\\User\\Documents\\circos-0.67-5\\test'

	directory = 'test_2'
	destination_directory = 'test_2'
	image_direct = 'test_2_circos'

	karyotype_name = 'karyotype.txt'
	link_output = 'links.txt'
	# TODO: replace with call to progrssiveMauve and taking test2_0409.fa from
	# the command line
	backbone_filename = 'test2_0409.xmfa.backbone'
	seq_filename = 'test2_0409.fa'
	output_conf_filename = 'test2_0409.conf'
	xmfa = os.path.join(directory, 'test2_0409.xmfa')
	backbone_file = os.path.join(directory, backbone_filename)
	alignment_list = xmfa_parse(xmfa)

	# Print alignments
	#for i in range(0, len(alignment_list[1]), 100):
		#print "\n%s..%s" % (i, i + 100)
		#print '\n'.join(x[i:i + 100] for x in alignment_list[1:])
	link_dict = links(backbone_file, link_output)
	link_dict = add_pct_identity(link_dict, os.path.join(directory, seq_filename), alignment_list)
	print link_dict
	import sys
	sys.exit()
	write_link_file(karyotype(os.path.join(directory, seq_filename), karyotype_name, directory), link_dict, link_output, directory)
	sample_conf = open('sample_conf.conf', 'r')
	output_conf = open(os.path.join(directory, output_conf_filename), 'w')
	for line in sample_conf.readlines():
		i = 0
		while i!= len(line):
			if line[0] == '<':
				output_conf.write(str(line).strip()+'\n')
				break
			elif line[i] != ' ':
				i+=1
			else:
				if line[0:i] == 'karyotype':
					output_conf.write('karyotype = '+destination_directory + str(karyotype_name)+'\n')
				elif line[0:i] == 'file':
					output_conf.write('file = '+destination_directory+str(link_output)+'\n')
				elif line[0:i] == 'dir*':
					output_conf.write('dir* = '+image_direct+'\n')
				elif line[0:i] == 'color':
					break
				else:
					if line[0] != '#':
						output_conf.write(str(line).strip()+'\n')
				break
	sample_conf.close()
	output_conf.close()


if __name__ == '__main__':
	main()
