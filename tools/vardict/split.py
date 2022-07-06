import sys


fai = sys.argv[1]
chunk_size = int(sys.argv[2])
overlap = int(sys.argv[3])  # Base pairs
with open(fai, 'r') as infile:
    for line in infile:
        name = line.split('\t')[0]
        stop = int(line.split('\t')[1])
        start = 1
        while start < stop:
            start = max(1, start - overlap)
            print('\t'.join([name, str(start),
                             str(min(start + chunk_size, stop))]))
            start += chunk_size
