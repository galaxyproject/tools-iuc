import random
import os.path
import copy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


class Nucleic_Acid(object):
	class DNA:
		def __init__(self, length=1000, seq=None):
			"""DNA class takes a length parameter and generates a sequence of
			that length.

			Optionally it can take a sequence string from a previously
			generated sequence.
			"""
			self.list_of_changes = []
			self.nucl_seq = Seq(self.generate_fragment(length), generic_dna)

			if seq is not None:
				self.nucl_seq = Seq(seq, generic_dna)

		@classmethod
		def generate_fragment(cls, length):
			letters = ('A', 'C', 'T', 'G')
			return ''.join([random.choice(letters) for i in range(length)])

		def insertion(self, length, start):
			insertion_fragment = self.generate_fragment(length)
			self.nucl_seq = self.nucl_seq[0:start] + insertion_fragment + self.nucl_seq[start:]
			self.list_of_changes.append('i%s@%s' % (length, start))

		def deletion(self, length, start):
			self.nucl_seq = self.nucl_seq[0:start] + self.nucl_seq[start + length:]
			self.list_of_changes.append('d%s@%s' % (length, start))

		def inversion(self, length, start):
			self.nucl_seq = self.nucl_seq[0:start] + \
				self.nucl_seq[start:start + length].reverse_complement() + \
				self.nucl_seq[start + length:]
			self.list_of_changes.append('v%s@%s' % (length, start))

		def translocation(self, length, start):
			translocated_segment = self.nucl_seq[start:start + length]
			new_DNA = self.nucl_seq[0:start] + self.nucl_seq[start + length:]
			new_start_pos = random.randint(0, len(new_DNA))
			self.nucl_seq = new_DNA[0:new_start_pos] + \
				translocated_segment + \
				new_DNA[new_start_pos:]
			self.list_of_changes.append('t%s@%s2%s' % (length, start, new_start_pos))

		def mutation_event(self):
			length = random.randint(1, 100)
			start = random.randint(0, len(self.nucl_seq))
			print length, start

			# List of methods to possibly call
			methods = ('insertion', 'deletion', 'translocation', 'inversion')
			# Randomly choose one of the above, and get that attribute of this Class.
			# getattr(self) lets you pick out a function by name and store in method_to_call
			method_to_call = getattr(self, random.choice(methods))
			# Which we then call and pass apporpriate parameters
			method_to_call(length, start)

if __name__ == '__main__':
	with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "test2_0409.fa"), 'w') as output:
		parent = Nucleic_Acid.DNA(random.randint(1000, 2000))

		parent_record = SeqRecord(parent.nucl_seq, id="Parent", description="length=%s" % len(parent.nucl_seq))
		SeqIO.write(parent_record, output, 'fasta')

		for i in range(5):
			obj2 = Nucleic_Acid.DNA(seq=str(parent.nucl_seq))
			# Use copy.deepcopy in order to copy the entire object. Simply
			# assigning obj2.nucl_seq = parent.nucl_seq is likely to have bad
			# effects.
			#
			# However, instead of copying, we'll just add an optional parameter
			# that will take a sequence in, avoiding any copying problems.
			for event in range(2):
				obj2.mutation_event()
			# Write obj2 to the fasta file
			SeqIO.write(SeqRecord(obj2.nucl_seq, id="Mut%s" % i, description=" ".join(obj2.list_of_changes)), output, 'fasta')
