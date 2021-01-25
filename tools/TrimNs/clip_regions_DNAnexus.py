import os
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Alphabet.IUPAC import ambiguous_dna
import gzip
import sys

def main():

	fasta_file = sys.argv[1]
	contamination_file = sys.argv[2]
	decontaminated_fasta_file = sys.argv[3]
	
	test_contamination_parsing = False
	verbose = True

	coord_list_for_sequence = {}

	mito_flag = False
	common_euk_flag = False

	# MASK always converts to Ns
	# TRIM always trims (trim_Ns usually results in this)
	# CONTAMINANT converts to Ns unless it is at the very start or end, in which case it is trimmed (Vecscreen is CONTAMINANT)
	# MITOCHONDRIAL is like CONTAMINANT, but the sequence is set aside in a special mito file
	# REMOVE removes
	# Note that a REMOVE in the mitochondrial section results in both a REMOVE and a MITOCHONDRIAL entry

	with open(contamination_file, 'r') as coord_input_handle:
		for line in coord_input_handle:
			if re.search('^\S*(REMOVE|TRIM|MASK|CLIP|VecScreen)', line):
				fields = re.split('\s+', line)
				id = fields[1]

				if re.search('^REMOVE', fields[0]):
					treatment = 'REMOVE'
					if id not in coord_list_for_sequence:
						coord_list_for_sequence[id] = []
					if mito_flag: # If this is mito sequence, flag it for copying too
						coord_list_for_sequence[id].append([0,0, 'MITOCHONDRIAL_REMOVE'])
					else:
						coord_list_for_sequence[id].append([0,0, treatment])

				else:
					start = fields[2]
					end = fields[3]

					treatment = 'MASK'
					if re.search('^(TRIM|REVCLIP|FWDCLIP|CLIP)', fields[0]):
						treatment = 'TRIM'
					if re.search('^VecScreen', fields[0]):
						treatment = 'CONTAMINANT'
					
					if id not in coord_list_for_sequence:
						coord_list_for_sequence[id] = []
					coord_list_for_sequence[id].append([int(start), int(end), treatment])
			elif (mito_flag or common_euk_flag) and not re.search('^#', line) and not re.search('^\=', line):
				fields = re.split('\s+', line)

				if len(fields) > 6:
					
					treatment = 'CONTAMINANT'
					if mito_flag:
						treatment = 'MITOCHONDRIAL' # CAPTURE MITO
					id = fields[0]
					start = int(fields[6])
					end = int(fields[7])
					start, end = sorted([start, end])
					if id not in coord_list_for_sequence:
						coord_list_for_sequence[id] = []
					coord_list_for_sequence[id].append([start, end, treatment])
			elif re.search('\=\=\=', line): # If this is a new section, this must be the end of the mito region
				mito_flag = False
				common_euk_flag = False
			if re.search('MITO', line): # If this is the mito region, change behaviour
				mito_flag = True
			if re.search('COMMON CONTAMINANTS IN EUKARYOTES', line): # If this is the mito region, change behaviour
				common_euk_flag = True

	# Merge entries
	#merged_coord_list_for_sequence = {}
	for id in coord_list_for_sequence:
		max_end = -100
		max_treatment = None
		for coord_pair_and_treatment in sorted(coord_list_for_sequence[id], key = lambda cpt: cpt[0]):
			if coord_pair_and_treatment[0] <= max_end:
				print('Screening coordinate overlap in ' + id + ' TREATMENT:' + coord_pair_and_treatment[2] + ' vs ' + max_treatment)
			if coord_pair_and_treatment[1] > max_end:
				max_end = coord_pair_and_treatment[1]
				max_treatment = coord_pair_and_treatment[2]

		termini = []
		coord_pair_and_treatment_for_label = {}

		for coord_pair_and_treatment in sorted(coord_list_for_sequence[id], key = lambda cpt: cpt[0]):
		
			# Reject any malformed features
			if coord_pair_and_treatment[0] > coord_pair_and_treatment[1] or coord_pair_and_treatment[0] < 0 or coord_pair_and_treatment[1] < 0:
				print('Malformed feature:', str(coord_pair_and_treatment[0]), '-', str(coord_pair_and_treatment[1]))
				next

			label = coord_pair_and_treatment[2] + ':' + str(coord_pair_and_treatment[0]) + '-' + str(coord_pair_and_treatment[1])
			coord_pair_and_treatment_for_label[label] = coord_pair_and_treatment

			start_position = {
				'TERMINUS': 'START',
				'POSITION': coord_pair_and_treatment[0],
				'LABEL': label,
			}

			end_position = {
				'TERMINUS': 'END',
				'POSITION': coord_pair_and_treatment[1],
				'LABEL': label,
			}			

			termini.append(start_position)
			termini.append(end_position)
			
		termini = sorted(termini, key=sort_termini)

		depth = 0
		current_overlap_labels = []
		overlap_groups = []

		for terminus in termini:
			if terminus['TERMINUS'] == 'START':
				depth += 1
				current_overlap_labels.append(terminus['LABEL'])
			else:
				depth -= 1
				if depth == 0:
					overlap_groups.append(current_overlap_labels)
					current_overlap_labels = []


		merged_coord_pairs_and_treatments = []
		merge_flag = False

		for overlap_group in overlap_groups:

			# If this label doesn't overlap, pass it to the merged set
			if len(overlap_group) == 1:
				merged_coord_pairs_and_treatments.append(coord_pair_and_treatment_for_label[overlap_group[0]])
			else:
				merge_flag = True
				starts = []
				ends = []
				treatments = {}
				for overlapping_label in overlap_group:
					#print(overlapping_label, end=' ')
					starts.append(coord_pair_and_treatment_for_label[overlapping_label][0])
					ends.append(coord_pair_and_treatment_for_label[overlapping_label][1])
					treatments[coord_pair_and_treatment_for_label[overlapping_label][2]] = 1

				starts = sorted(starts)
				ends = sorted(ends)
				min_start = starts[0]
				max_end = ends[-1]

				if len(treatments) == 1:
					merged_coord_pair_and_treatment = [
						min_start,
						max_end,
						list(treatments.keys())[0],
					]
					merged_coord_pairs_and_treatments.append(merged_coord_pair_and_treatment)
				elif len(treatments) ==2 and 'TRIM' in treatments and 'CONTAMINANT' in treatments:
					merged_coord_pair_and_treatment = [
						min_start,
						max_end,
						'TRIM',
					]
					merged_coord_pairs_and_treatments.append(merged_coord_pair_and_treatment)
				else:
					quit('CANNOT MERGE: ' + str(overlap_group))

		if merge_flag:
			print('UNMERGED: ' + str(coord_list_for_sequence[id]))
			coord_list_for_sequence[id] = merged_coord_pairs_and_treatments
			print('  MERGED: ' + str(merged_coord_pairs_and_treatments))

	if test_contamination_parsing:
		quit('Test run')

	mask_flag = False

	fasta_output_handle = open(decontaminated_fasta_file, 'w')

	mitochondrial_records = []

	with open(fasta_file, "r") as fasta_input_handle:
		for record in SeqIO.parse(fasta_input_handle, "fasta"):

			remove_flag = False

			if(record.id in coord_list_for_sequence):
				if verbose:
					print('Extracting from', record.id)

				last_end = 0
				edited_record = SeqRecord(Seq('', ambiguous_dna), id=record.id, description = record.description, name = record.name)
				for coord_pair_and_treatment in sorted(coord_list_for_sequence[record.id], key = lambda cpt: cpt[0]):
	
					if verbose:
						print('\t', coord_pair_and_treatment)

					# We handle removal of mito sequence by treating it as a REMOVE entry and a MITOCHONDRIAL entry
					# The MITOCHONDRIAL entry initially has coords -1,-1 until this stage, where we know the length
					if coord_pair_and_treatment[2] == 'MITOCHONDRIAL_REMOVE' and coord_pair_and_treatment[1] < 1:
						coord_pair_and_treatment[0] = 1
						coord_pair_and_treatment[1] = len(record.seq)

					if coord_pair_and_treatment[2] == 'REMOVE' or coord_pair_and_treatment[2] == 'MITOCHONDRIAL_REMOVE':
						remove_flag = True

					if coord_pair_and_treatment[2] != 'REMOVE':
						edited_record.seq += record.seq[last_end:(coord_pair_and_treatment[0]-1)] # Takes account of 0-based numbering
						if coord_pair_and_treatment[2] == 'MASK' or ((coord_pair_and_treatment[2] == 'CONTAMINANT' or coord_pair_and_treatment[2] == 'MITOCHONDRIAL') and (coord_pair_and_treatment[0] > 1 and coord_pair_and_treatment[1] < len(record.seq))):
							edited_record.seq += 'N' * (coord_pair_and_treatment[1] - coord_pair_and_treatment[0] + 1)
						if coord_pair_and_treatment[2] == 'MITOCHONDRIAL' or coord_pair_and_treatment[2] == 'MITOCHONDRIAL_REMOVE':
							mitochondrial_record = SeqRecord(Seq('', ambiguous_dna), id=str(len(mitochondrial_records)+1),description='', name = '')
							mitochondrial_record.seq += record.seq[(coord_pair_and_treatment[0]-1):(coord_pair_and_treatment[1])]
							mitochondrial_records.append(mitochondrial_record)

						last_end = coord_pair_and_treatment[1]

				edited_record.seq += record.seq[last_end:] # Finally add the end segment
				record = edited_record

			if not(remove_flag):
				SeqIO.write([record], fasta_output_handle, 'fasta')

	fasta_output_handle.close()

def sort_termini(terminus):
	if terminus['TERMINUS'] == 'END':
		terminus['POSITION'] += 0.5
	return terminus['POSITION']

if __name__ == "__main__":
		main()
