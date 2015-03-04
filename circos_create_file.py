def links(backbone_file,link_output):
	link_output = open(link_output,'w')
	backbone_file = open(backbone_file,'r')
	new_readlines = []
	for line in backbone_file.readlines():
	    t_ct = 0
	    new_line = 'chr1 '
	    for character in line:
	        if character == '\t':
	            t_ct+=1
	            if t_ct == 2:
	                new_line += ' chr2 '
	            elif t_ct == 4:
	                new_line += ' chr3 '
	            else:
	                new_line += ' '
	        elif character == '\n':
	            pass
	        else:
	            new_line += character
	    if new_line.count('0 0')<>2:
	        new_readlines += [new_line]
	    else:
	        pass
	for line in new_readlines:
	   new_line_1 = ''
	   new_line_2 = ''
	   new_line_3 = ''
	   line = line.replace('chr1 0 0',"")
	   line = line.replace('chr2 0 0',"")
	   line = line.replace('chr3 0 0',"")
	   line = line.strip()
	   if 'chr1' in line and 'chr2' in line and 'chr3' in line:
	        chr2_pos = line.find('chr2')
	        chr3_pos = line.find('chr3')
	        new_line_1 = (line[chr2_pos-1:]).strip()
	        new_line_2 = (line[:chr3_pos-1]).strip()
	        new_line_3 = (line[:chr2_pos-1]+line[chr3_pos-1:]).strip()
	   if line.count('leftend')==0 and new_line_1 == '':
	        link_output.write('%s\n' %line)
	   elif line.count('leftend')==0:
	        link_output.write(new_line_1+str('\n'))
	        link_output.write(new_line_2+str('\n'))
	        link_output.write(new_line_3+str('\n'))
	backbone_file.close()
	link_output.close()


def karyotype(seq_file):
	genome_list = []
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
			name = line.strip()
	seq_file.close()
	return genome_list

def main():
	backbone_file = ''
	link_output = ''
	sample_conf = open('sample_conf.conf','r')
        output_conf = open('output.conf','w')
	for line in sample_conf.readlines():
			output_conf.write(str(line).strip()+'\n')
	sample_conf.close()
	genome_list = karyotype('second_run.fa')
	print genome_list


	output_conf.close()

if __name__ == '__main__':
	main()

