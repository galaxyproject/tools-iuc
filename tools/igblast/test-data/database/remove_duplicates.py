#!/usr/bin/env python


'''
	I have seen that in one working copy of the database 
	(/mnt/home/aerijman/lib/miniconda3/envs/presto6/share/igblast/database-6d21f2f0-9860-4b79-a9e9-193f805fbf73/v.fa),
	the fastas are all UPPER CASE and NOT WRAPED (every 60 characters), hence, I add here such features (at the very end). 
'''


import sys
import numpy as np

# get filename from argument
for n, i in enumerate(sys.argv):
    if i in ['-in','--in','-input','--input']:
        fasta = sys.argv[n+1]

if not 'fasta' in locals():
    sys.exit('python <this script> -in <file.fasta>')

data, seq = {}, ""

# check if repeated key corresponds to repeated value
check_status = False

# open file
with open(fasta) as f:
    while True:
        try:
            line = next(f).strip()

            if len(line)<2:
                continue            
   
            # case header
            if line[0] == ">":
                if len(seq)>0:  # first line doesn't have seq yet
    
                    # don't repeat data
                    if header in data:  
                        if seq != data[header]:
                            sys.stderr.write("found repeated header with different sequences: {}\n".format(header))
                        else:
                            sys.stderr.write("avoiding re-writing same data more than once for {}\n".format(header))
                    else:           
                        data[header] = seq
                
                seq = ""
                header = line[1:]  # esclude ">" symbol


            # case sequence
            else:
                seq += line
                
        except StopIteration:

            # don't repeat data
            if check_status and seq==data[header]:
                print("found repeated sequence: {}".format(header))
            else:
                data[header] = seq

            break

for k,v in data.items():
    v_fasta = v.upper()  #'\n'.join([v[i:i+60] for i in np.arange(0, len(v), 60)])
    print(">{}\n{}".format(k, v_fasta))
