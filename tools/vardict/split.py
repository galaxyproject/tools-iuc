import sys


fai = sys.argv[1]
chunk_size = int(sys.argv[2])
overlap = 150 # Base pairs
with open(fai, 'r') as infile:
    for line in infile:
        name = line.split('\t')[0]
        stop = int(line.split('\t')[1])
        start = 1
        while start < stop:
            start = max(1, start-overlap)
            out = '\t'.join([name, str(start), str(min(start+chunk_size, stop))])
            print(out)
            # outfile.write(out + '\n')
            start += chunk_size
