import random
def insertion(mut_seq,changes_list):
	i = random.randint(0,len(mut_seq)-1)
	insertion_length = int(random.uniform(1,100))
	#Probability should decrease with insertion size?
	for number in range(insertion_length):
		bases = ('A','T','C','G')
		base = random.choice(bases)
		new_seq = mut_seq[:i] + str(base) + mut_seq[i:]
		changes_list += ['Insertion of '+str(base)+' at position '+str(i+1)]
	return new_seq,changes_list

def rev_comp(sequence):
	pairing_dict = {'A':'T','G':'C','T':'A','C':'G'}
	nucl_list = [pairing_dict[letter] for letter in sequence[::-1]]
	return ''.join(str(char) for char in nucl_list)

def deletion(mut_seq,changes_list):
	i = random.randint(0,len(mut_seq)-1)
	j = random.randint(i+1,i+100)
	if j>=len(mut_seq):
		j = len(mut_seq)
	new_seq = mut_seq[:min(i,j)] + mut_seq[max(i,j):]
	changes_list += ['Deletion at position '+str(i+1)]
	return new_seq,changes_list

def inversion(mut_seq,changes_list):
	i = random.randint(0,len(mut_seq)-1)
	j = random.randint(0,len(mut_seq)-1)
	if i==j:
		return mut_seq,changes_list
	elif i<j:
		endpoints = (i,j)
	else:
		endpoints = (j,i)
	new_seq = mut_seq[:i]+rev_comp(mut_seq[i:j])+mut_seq[j:]
	changes_list += ['Inversion between positions ' + str(endpoints[0]+1) + ' and ' + str(endpoints[1])]
	return new_seq,changes_list

def translocation(mut_seq,changes_list):
	i = random.randint(0,len(mut_seq)-1)
	j = random.randint(0,len(mut_seq)-1)
	if i==j:
		changes_list += 'No Change'
		return mut_seq,changes_list
	elif i<j:
		endpoints = (i,j)
	else:
		endpoints = (j,i)
	transposon = mut_seq[endpoints[0]:endpoints[1]]
	excise_seq = mut_seq[0:endpoints[0]]+mut_seq[endpoints[1]:]
	k = random.randint(0,len(excise_seq)-1)
	l = k+len(transposon)
	if k==l:
		1==1
	elif k<l:
		new_pos = (k,l)
	else:
		new_pos = (l,k)
	new_seq = excise_seq[0:new_pos[0]]+transposon+excise_seq[new_pos[1]-len(transposon):]
	changes_list += ['Translocation from (' + str(endpoints[0]+1) + ',' + str(endpoints[1]) + ') to ' +str(new_pos)]
	return new_seq,changes_list




def main():
	number_of_sequences = 2
	change_list_output = open('3_genome_data/three_genomes.txt', 'w')
	sequence_output = open('3_genome_data/three_genomes.fa', 'w')
	bases = ('A','T','C','G')
	mutation_probs = [45,35,15,5]
	base_seq= ''
	for iteration in range(10000):
		base_seq += random.choice(bases)
	mut_seq = base_seq
	sequence_output.write('>original\n%s\n' % base_seq)
	i=1
	for sequence in range(number_of_sequences):
		changes_list = []
		for evolutionary_event in range(500):
			Prob = random.randint(0,100)
			if Prob in range(0,mutation_probs[0]):
				mut_seq,changes_list = insertion(mut_seq,changes_list)
			elif Prob in range(mutation_probs[0],mutation_probs[0]+mutation_probs[1]):
				mut_seq,changes_list = deletion(mut_seq,changes_list)
			elif Prob in range(sum(mutation_probs[0:2]),sum(mutation_probs[0:3])):
				mut_seq,changes_list = inversion(mut_seq,changes_list)
			else:
				mut_seq,changes_list = translocation(mut_seq,changes_list)
		sequence_output.write('>mut'+str(i)+'\n%s\n' % mut_seq)
		for change in changes_list:
			change_list_output.write(change + "\n")
		change_list_output.write("End of Sequence" + str(i) + '\n')
		sequence_output.close
		i+=1
	change_list_output.close()


if __name__ == '__main__':
	main()
