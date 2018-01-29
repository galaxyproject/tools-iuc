from __future__ import print_function

import argparse
import fileinput
import sys

# Maximum value of a signed 32 bit integer (2**31 - 1).
MAX_CHROM_LEN = 2147483647


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)


parser = argparse.ArgumentParser()
parser.add_argument('--input', dest='input', help="Input dataset")
parser.add_argument('--start_coordinate', dest='start_coordinate', type=int, help='Chromosome start coordinate, either 0 or 1.')
parser.add_argument('--subtract_from_start', dest='subtract_from_start', type=int, help='Distance to subtract from start.')
parser.add_argument('--add_to_end', dest='add_to_end', type=int, help='Distance to add to end.')
parser.add_argument('--extend_existing', dest='extend_existing', help='Extend existing start/end instead of from computed midpoint.')
parser.add_argument('--chrom_len_file', dest='chrom_len_file', help="File names of .len files for chromosome lengths")
parser.add_argument('--region_boundaries', dest='region_boundaries', help="Option for handling region boundaries")
parser.add_argument('--output', dest='output', help="Output dataset")
args = parser.parse_args()

extend_existing = args.extend_existing == 'existing'
out = open(args.output, 'wb')

chrom_start = int(args.start_coordinate)
chrom_lens = dict()
# Determine the length of each chromosome and add it to the chrom_lens dictionary.
len_file_missing = False
len_file_error = None
len_file = fileinput.FileInput(args.chrom_len_file)
try:
    for line in len_file:
        fields = line.split("\t")
        chrom_lens[fields[0]] = int(fields[1])
except Exception as e:
    len_file_error = str(e)

with open(args.input) as fhi:
    for line in fhi:
        if line.startswith('#'):
            # Skip comments.
            continue
        items = line.split('\t')
        if len(items) != 9:
            # Skip invalid gff data.
            continue
        chrom = items[0]
        start = int(items[3])
        end = int(items[4])
        if extend_existing:
            new_start = start - args.subtract_from_start
            new_end = end + args.add_to_end
        else:
            midpoint = (start + end) // 2
            new_start = midpoint - args.subtract_from_start
            new_end = midpoint + args.add_to_end
        # Check start boundary.
        if new_start < chrom_start:
            if args.region_boundaries == 'discard':
                continue
            elif args.region_boundaries == 'limit':
                new_start = chrom_start
            elif args.region_boundaries == 'error':
                out.close()
                stop_err('Requested expansion places region beyond chromosome start boundary of %d.' % chrom_start)
        # Check end boundary.
        chrom_len = chrom_lens.get(chrom, None)
        if chrom_len is None:
            len_file_missing = True
            chrom_len = MAX_CHROM_LEN
        if new_end > chrom_len:
            if args.region_boundaries == 'discard':
                continue
            elif args.region_boundaries == 'limit':
                new_end = chrom_len
            elif args.region_boundaries == 'error':
                out.close()
                stop_err('Requested expansion places region beyond chromosome end boundary of %d.' % chrom_len)
        new_line = '\t'.join([chrom, items[1], items[2], str(new_start), str(new_end), items[5], items[6], items[7], items[8]])
        out.write(new_line)
out.close()

if len_file_error is not None:
    print("All chrom lengths set to %d, error in chrom len file: %s" % (MAX_CHROM_LEN, len_file_error))
if len_file_missing:
    print("All chrom lengths set to %d, chrom len files are not installed." % MAX_CHROM_LEN)
