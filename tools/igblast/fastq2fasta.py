#!/usr/bin/env python

import sys
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

if '-in' not in sys.argv:
	sys.exit('\n python <this script> -in <fastq file>\n')

for n,i in enumerate(sys.argv):
	if i=='-in':
		fastq = sys.argv[n+1]

fasta=[]

with open(fastq, 'r') as f:
	while True:
		try:
			a,b,c,d = next(f), next(f), next(f), next(f)
			fasta.append('>{}'.format(a) + b)

		except StopIteration:
			break
		
		n+=1


print(''.join(fasta))
