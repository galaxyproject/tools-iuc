import os
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Alphabet.IUPAC import ambiguous_dna
import sys

def main():

	fasta_file = sys.argv[1]
	output_file = sys.argv[2]

	minleftover = 100   # after trimming start/end, at least this many bp should be left
	winsize     = 5000  # for sliding window analysis
	minslidingBase = 0.4   # maximum fraction of Ns in sliding window before alarm sets off

	output_handle = open(output_file, 'w')
	with open(fasta_file, "r") as fasta_input_handle:
		for record in SeqIO.parse(fasta_input_handle, "fasta"):
			# output_handle.write('Testing ' + record.id + "\n")

			# Do this twice: once with the regular string, once with the reversed string
			# The output should be as one or more variables that can be reversed for the second iteration
			seq_string = str(record.seq)

			n_count = seq_string.count('N') + seq_string.count('n')
			n_perc = n_count/ len(seq_string)
			if n_perc > 0.8:
				output_handle.write('# WARNING: ' + record.id + '\t' + str(int(n_perc * 10000)/100) + ' % Ns of total ' + str(len(seq_string)) + "\n")

			start_n_match = re.match('^([Nn]*)', seq_string)
			end_n_match = re.search('^([Nn]*)', seq_string[::-1])

			startseq = ''
			if start_n_match:
				startseq = start_n_match.group(1)
			endseq = ''
			if end_n_match:
				endseq = end_n_match.group(1)
			realseq_length = len(seq_string) - len(startseq) - len(endseq)
			# Handle "all Ns exception"
			if len(startseq) == len(seq_string):
				realseq_length = 0
				endseq = '' 

			if len(startseq) > 0 or len(endseq) > 0:
				if len(startseq) > 0:
					output_handle.write('\t'.join(["TRIM:",record.id,'1',str(len(startseq))]) + "\t\n")
				if len(endseq) > 0:
					output_handle.write('\t'.join(["TRIM:",record.id,str(len(startseq)+realseq_length+1),str(len(seq_string))]) + "\t\n")

			if realseq_length <= minleftover:
				output_handle.write("REMOVE: " + record.id + '\t' + str(realseq_length) + " bp leftover after trimming\n")
			
			# Only attempt the windowing approach if there's an initial window that doesn't meet the trigger condition
			for seq_string_loop, seq_string_for_window in enumerate([seq_string, seq_string[::-1]]):
				if len(seq_string_for_window) > winsize and (seq_string_for_window[:winsize].count('N') + seq_string_for_window[:winsize].count('n')) > winsize * minslidingBase:

					non_n_regions = []
					non_n_iterator = re.finditer('[^Nn]+', seq_string_for_window)
					for non_n_instance in non_n_iterator:

						current_non_n_start = non_n_instance.start(0) + 1
						current_non_n_end = non_n_instance.end(0)

						non_n_regions.insert(0, [current_non_n_start, current_non_n_end])

						# Does the *end* of this block satisfy the window condition?
						bases_in_window = 0
						start_of_first_base_block_in_window = None
						for non_n_region in non_n_regions:
							if non_n_region[1] >= current_non_n_end - winsize:
								start_of_first_base_block_in_window = non_n_region[0]
								if non_n_region[0] >= current_non_n_end - winsize:
									bases_in_window += non_n_region[1] - non_n_region[0] + 1
								else:
									bases_in_window += non_n_region[1] - non_n_region[0] + 1
									break
							else:
								break

						if bases_in_window >= minslidingBase * winsize:

							# Remember that the blocks are in *reverse* order along the sequence
							trackback_to_start_flag = False

							for i, non_n_region in enumerate(non_n_regions):
								if i == 0:
									continue
								bases_in_window_2 = 0
								if non_n_region[1] < non_n_regions[0][0] - winsize:
									break
								else:
									current_window_start = max(non_n_region[0], non_n_regions[0][0] - winsize)

									# Count the bases from this block
									bases_in_window_2 += min(non_n_region[1], (current_window_start+winsize-1) ) - current_window_start + 1

									# Add up all the sequence in blocks tested thus far, but not the final block
									for j in range(1, i):
										bases_in_window_2 += non_n_regions[j][1] - non_n_regions[j][0] + 1

									# Count the sequence in the final block that can contribute
									# This is the region between the start of the final block and either the end of the block
									# or the end of a window extending from the current test start point, whichever comes first
									bases_in_window_2 += min( non_n_regions[0][1], (current_window_start + winsize-1) ) - non_n_regions[0][0] + 1

									# output_handle.write('BIW: ' + str(i) + ' ' + str(bases_in_window_2) + "\n")
									if bases_in_window_2 >= minslidingBase * winsize:
										if current_window_start == non_n_region[0]:
											start_of_first_base_block_in_window = current_window_start
										else:
											start_of_first_base_block_in_window = non_n_regions[i-1][0]
									else:
										# We keep going back until the threshold condition is *not* met
										break

									# If we find the break-point should be before the first block, then we don't want to trim at all!										
									if i == len(non_n_regions) - 1:
										trackback_to_start_flag = True


							# Only trim if the breakpoint isn't before the first block
							# and if the breakpoint isn't at the start of the sequence
							if not(trackback_to_start_flag) and start_of_first_base_block_in_window != 1:
								if seq_string_loop == 0:
									output_handle.write("FWDCLIP:\t" + record.id + "\t1\t" + str(start_of_first_base_block_in_window - 1) + "\n")
								else:
									output_handle.write("REVCLIP:\t" + record.id + '\t' + str(len(seq_string) - start_of_first_base_block_in_window + 2) + '\t' + str(len(seq_string)) + "\n")

							break

if __name__ == "__main__":
		main()
