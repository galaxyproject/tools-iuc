import pprint
import os.path
import itertools
import re
import subprocess
import tempfile
import argparse
from Bio import SeqIO


def xmfa_parse(xmfa_handle):
	start_regex = re.compile('> (?P<idx>[1-9]+):(?P<start>[0-9]+)-(?P<end>[0-9]+)')
	seq_info_list = {}
	seqidx = None
	seq_start = None

	for line in xmfa_handle.readlines():
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


def links(backbone_file_handle):
	header_parsed = False
	links = {}
	for line in backbone_file_handle:
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


def write_link_file(names_list, links, link_output):
	with open(link_output,  'w') as handle:
		for key in links:
			key_from, key_to = key.split('-')  # Rememver, we keyed on 1-3 0-2 3-2 / etc
			for link in links[key]:  # List of from/to
				print link
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



def karyotype(seq_file):
	# TODO: brewer colours
	colors_list = ['red', 'blue', 'green', 'orange', 'violet', 'brown']
	genome_list = []
	karyotype = []

	for i, record in enumerate(SeqIO.parse(seq_file, 'fasta')):
		genome_list.append(record.id)
		karyotype.append('chr - %s %s 0 %s %s\n' %
						 (record.id, record.id, len(record.seq), colors_list[i]))
	return genome_list, karyotype


def add_pct_identity(link_dict, alignment_list):
	for key in link_dict:
		key_pair = map(int, key.split('-'))

		good_links = []
		for element in link_dict[key]:
			first_index_set = [
				min(abs(element[0]), abs(element[1])),
				max(abs(element[0]), abs(element[1]))]
			second_index_set = [
				min(abs(element[2]), abs(element[3])),
				max(abs(element[2]), abs(element[3]))]

			sect_1 = alignment_list[key_pair[0]][first_index_set[0]:first_index_set[1]]
			sect_2 = alignment_list[key_pair[1]][second_index_set[0]:second_index_set[1]]

			if element[0] < 0:
				sect_1 = reverse_complement(sect_1)
			if element[2] < 0:
				sect_2 = reverse_complement(sect_2)

			if first_index_set[0] != 0 and second_index_set != 0:
				good_links.append(element + [percent_sequence_identity(sect_1, sect_2)])

		link_dict[key] = good_links
	pprint.pprint(link_dict)
	return link_dict


def main():
	parser = argparse.ArgumentParser(description='Circos conf from progressiveMauve output')
	parser.add_argument('xmfa', type=file, help='ProgressiveMauve XMFA file')
	parser.add_argument('backbone', type=file, help='ProgressiveMauve backbone file')
	parser.add_argument('sequence', type=file, help='Fasta sequence')
	parser.add_argument('--output', help="Output Directory", default="test_2")
	parser.add_argument('--verbose', action='store_true', help="Verbose output/logging")
	args = parser.parse_args()

	# Create it if it doesn't exist
	if not os.path.exists(args.output):
		os.mkdir(args.output)

	# If not a directory, die early
	if not os.path.isdir(args.output):
		raise Exception("--output must be a directory")

	destination_directory = 'C:\\Users\\User\\Desktop\\491 Scripts\\test_2\\'
	image_direct = 'C:\\Users\\User\\Documents\\circos-0.67-5\\test'

	destination_directory = 'test_2'
	image_direct = 'test_2_circos'

	# TODO: replace with call to progrssiveMauve and taking test2_0409.fa from
	# the command line
	output_conf_filename = 'test2_0409.conf'

	alignment_list = xmfa_parse(args.xmfa)

	# Print alignments
	if args.verbose:
		for i in range(0, len(alignment_list[1]), 100):
			print "\n%s..%s" % (i, i + 100)
			print '\n'.join(x[i:i + 100] for x in alignment_list[1:])

	# Construct genome-genome link dictionary
	link_dict = links(args.backbone)
	# Add %ID to each of the links
	link_dict = add_pct_identity(link_dict, alignment_list)
	# Collect 'chromosome' IDs and ...
	genome_id_list, karyotype_data = karyotype(args.sequence)
	# write karyotype data to output/karyotype.txt
	with open(os.path.join(args.output, 'karyotype.txt'), 'w') as handle:
		for chrom in karyotype_data:
			handle.write(chrom)
	# Write output/links.txt
	write_link_file(genome_id_list, link_dict,
					os.path.join(args.output, 'links.txt'))

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
