import pprint
import os.path
import itertools
def links(backbone_file,link_output):
	lines = open(backbone_file,'r').readlines()
	# There will be 2N where N is the number of genomes in our header
	header = lines[0].split('\t') 
	genome_count = len(header) / 2
	# Data structure to hold the links
	links = {}
	# A from-to list (from col/to col) which will be used to check links between genomes
	# https://docs.python.org/2/library/itertools.html
	# from_to for a genome_count of 3 looks like: [(0, 1), (0, 2), (1, 2)]
	# Preferred as it will expand to ANY number of comparisons
	from_to = list(itertools.combinations(range(genome_count), 2))
	# Now down to actual parsing
	for line in lines[1:]:
	  # Split with a tab. I can't remember if the file is tab separated, if not, then use `.split()` and it'll split on any whitespace.
	  link_data = line.split('\t')
	  # Iterate over pairs of columns
	  for x, y in from_to:
	    # Access the link data in that column
	    a_left = link_data[2*x].strip()
	    a_right = link_data[2*x + 1].strip()
	    b_left = link_data[2*y].strip()
	    b_right = link_data[2*y + 1].strip()
	    # if any of them are zero, then we can continue, as this isn't a "true" link.
	    # Circos will plot links from "0 0", so any links with "0 0" need to be removed.
	    if a_left != 0  and b_left != 0:
	        try:
	            links['%s-%s' % (x,y)].append([a_left, a_right, b_left, b_right])
	        except:
	            links['%s-%s' % (x,y)] = [
	          [a_left, a_right, b_left, b_right]
		    ]
	return links

def write_link_file(names_list,links, link_output='links.txt'):
  with open(os.path.join('/home/users/CPT/CPT/491/scrosby/Circos/3_genome_data',link_output), 'w') as handle:
    for key in links:
      key_from, key_to = key.split('-') # Rememver, we keyed on 1-3 0-2 3-2 / etc
      for link in links[key]: # List of from/to
        # Create the list by re-arranging the link_data
        # Circos will plot links from "0 0", so any links with "0 0" need to be removed.
	if link[0] != link[1] and link[2] != link[3]: 
            data = [names_list[int(key_from)]] + link[0:2] + [names_list[int(key_to)]] + link[2:4]
            handle.write(' '.join(data) + "\n")

def karyotype(seq_file):
	genome_list = []
	names_list = []
	seq_file = open(seq_file,'r')
	preceding_line_was_header = 0
	for line in seq_file.readlines():
                if preceding_line_was_header == 1:
			genome_tuple = (name,len(line))
			genome_list += [genome_tuple]
                preceding_line_was_header = 0
		name = ''
		genome_tuple = ('','')
                if line[0] == '>':
                        preceding_line_was_header = 1
			#Circos does not like ">" in the genome names.
			name = line[1:].strip()
	seq_file.close()
	karyotype_file = open(os.path.join('/home/users/CPT/CPT/491/scrosby/Circos/3_genome_data','3_genome_karyotype.txt'),'w')
	for tuple in genome_list:
		names_list+=[tuple[0]]
		karyotype_file.write('chr - '+(tuple[0]+' ')*2+str(0)+' '+str(tuple[1]-1)+' \n')
	karyotype_file.close()
	return names_list

def main():
	backbone_file = os.path.join('/home/users/CPT/CPT/491/scrosby/Circos/3_genome_data','3_genome_data.xmfa.backbone')
	link_output = '3_genome_links.txt'
	sample_conf = open('sample_conf.conf','r')
        output_conf = open(os.path.join('/home/users/CPT/CPT/491/scrosby/Circos/3_genome_data','three_genome_conf.conf'),'w')
	for line in sample_conf.readlines():
		if line[0] != '#':
			output_conf.write(str(line).strip()+'\n')
	sample_conf.close()
	output_conf.close()
	write_link_file(karyotype(os.path.join('/home/users/CPT/CPT/491/scrosby/Circos/3_genome_data','three_genome_data.fa')),links(backbone_file,link_output))
	
if __name__ == '__main__':
	main()

