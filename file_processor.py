link_output = open('link_coordinates.txt','w')
file = open('threeseq.xmfa.backbone','r')
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
   line = line.replace('chr1 0 0',"")
   line = line.replace('chr2 0 0',"")
   line = line.replace('chr3 0 0',"")
   line = line.strip()
   if line.count('leftend')==0:
        link_output.write('%s\n' %line)    
