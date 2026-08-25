import sys


fai = sys.argv[1]
chunk_size = int(sys.argv[2])
overlap = int(sys.argv[3])  # Base pairs
# vardictcpp requires a 4th (gene/name) column because -g 4 points past the
# coordinate columns; emit the contig name as the region label.
with open(fai, 'r') as infile:
    for line in infile:
        name = line.split('\t')[0]
        stop = int(line.split('\t')[1])
        start = 1
        while start < stop:
            start = max(1, start - overlap)
            print('\t'.join([name, str(start),
                             str(min(start + chunk_size, stop)), name]))
            start += chunk_size
