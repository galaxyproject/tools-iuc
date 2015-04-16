import random
import os.path

class Nucleic_Acid:
	class DNA:
		def __init__(self, length):
			self.variety = 'DNA'
			self.nucl_seq = ''
			self.list_of_changes = []
			for nucleotide in range(length):
				self.nucl_seq+=('A', 'T', 'G', 'C')[random.randint(0, 3)]

		def reverse_complement(self, start, end):
			DNA_pairing_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}
			rev_comp_template = self.nucl_seq[start:end]
			rev_comp_product = ''
			for letter in rev_comp_template[::-1]:
				if self.variety == 'DNA':
					rev_comp_product += DNA_pairing_dict[letter]
			return rev_comp_product

		def mutation_event(self):
			length = random.randint(1, 100)
			start = random.randint(0, len(self.nucl_seq))
			def insertion(self, length, start):
				insertion_fragment = ''.join([('A', 'T', 'G', 'C')[random.randint(0, 3)] for instance in range(length)])
				self.nucl_seq = self.nucl_seq[:start+1]+insertion_fragment+self.nucl_seq[start+1:]

			def deletion(self, length, start):
				self.nucl_seq = self.nucl_seq[:start+1]+self.nucl_seq[start+length+1:]

			def inversion(self, length, start):
				self.nucl_seq = self.nucl_seq[:start+1]+self.reverse_complement(start+1, start+1+length)+self.nucl_seq[start+1+length:]
				self.list_of_changes += ('inversion', start, start+length)
			def translocation(self, length, start):
				translocated_segment = self.nucl_seq[start:start+length+1]
				new_DNA = self.nucl_seq[:start]+self.nucl_seq[start+length+1:]
				new_start_pos = random.randint(0, len(new_DNA))
				self.nucl_seq = new_DNA[:new_start_pos]+translocated_segment+new_DNA[new_start_pos:]
			prob = random.randint(1, 100)
			if prob<=35:
				insertion(self, length, start)
			elif prob<=70:
				deletion(self, length, start)
			if prob<=90:
				translocation(self, length, start)
			else:
				inversion(self, length, start)



if __name__ == '__main__':
	output = open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "test2_0409.fa"), 'w')
	init_list = [Nucleic_Acid.DNA(random.randint(1000, 2000)) for number_of_sequences in range(1)]
	mut_list = []
	ops_list = []
	for number_of_sequences in range(5):
		obj2 = Nucleic_Acid.DNA(random.randint(1000, 2000))
		obj2.nucl_seq = init_list[0].nucl_seq
		for event in range(500):
			obj2.mutation_event()
		mut_list += [obj2]
		ops_list += [obj2.list_of_changes]
	output.write('>Original\n')
	output.write(str(init_list[0].nucl_seq)+'\n')
	i=0
	for sequence in mut_list:
		i+=1
		output.write('>Mut%s\n' % i)
		output.write(str(sequence.nucl_seq) + '\n')
	output.close()
   #print ops_list

