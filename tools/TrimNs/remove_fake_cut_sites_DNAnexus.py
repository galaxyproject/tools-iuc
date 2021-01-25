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
	log_file = sys.argv[3]

	output_handle = open(output_file, 'w')
	log_handle = open(log_file, 'w')

	with open(fasta_file, "r") as fasta_input_handle:
		for record in SeqIO.parse(fasta_input_handle, "fasta"):

			change_count = 0
			cut_sites = [
				Seq("CTTAAG"),
				Seq("CTTCTCG"),
				Seq("GCTCTTC"),
				Seq("CCTCAGC"),
				Seq("GAATGC"),
				Seq("GCAATG"),
				Seq("ATCGAT"),
				Seq("CACGAG"),
			]

			for cut_site in cut_sites:
				cut_site_both_orientations = ( cut_site, cut_site.reverse_complement() )

				for cut_site_for_orientation in cut_site_both_orientations:

					n_flank_length = 1
					search_pattern = 'N' * n_flank_length + str(cut_site_for_orientation) + 'N' * n_flank_length
					replacement = 'N' * (n_flank_length * 2 + len(cut_site_for_orientation))

					(new_string, changes) = re.subn(search_pattern, replacement, str(record.seq.upper()), flags=re.IGNORECASE)
					change_count += changes

					record.seq = Seq(new_string)

			if change_count > 0:
				log_handle.write(' '.join([record.id, ':', str(change_count), 'changes\n']))
			SeqIO.write([record], output_handle, 'fasta')

			# Finally, count the matches
			possible_fake_cut_sites = re.findall('N[^N]{1,10}N', str(record.seq.upper()))
			if len(possible_fake_cut_sites) > 0:
				log_handle.write(' '.join([record.id, ':', str(len(possible_fake_cut_sites)), 'possible non-standard fake cut sites\n']))

	output_handle.close()
	log_handle.close()

if __name__ == "__main__":
		main()
