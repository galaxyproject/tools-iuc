link_output = open('2nd_run_links.txt','w')
file = open('2nd_run.xmfa.backbone','r')
new_readlines = []
for line in file.readlines():
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
