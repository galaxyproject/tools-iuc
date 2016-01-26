import argparse
import sys


def stop_err( msg ):
    sys.stderr.write( msg )
    sys.exit(1)

parser = argparse.ArgumentParser()
parser.add_argument('--input', dest='input', help="Input dataset")
parser.add_argument('--subtract_from_start', dest='subtract_from_start', type=int, help='Distance to subtract from start.')
parser.add_argument('--add_to_end', dest='add_to_end', type=int, help='Distance to add to end.')
parser.add_argument('--extend_existing', dest='extend_existing', help='Extend existing start/end rather or from computed midpoint.')
parser.add_argument('--output', dest='output', help="Output dataset")
args = parser.parse_args()

extend_existing = args.extend_existing == 'existing'
out = open(args.output, 'wb')

for line in open(args.input):
    if line.startswith('#'):
        continue
    items = line.split('\t')
    if len(items) != 9:
        continue
    start = int(items[3])
    end = int(items[4])
    if extend_existing:
        start -= args.subtract_from_start
        end += args.add_to_end
    else:
        midpoint = (start + end) // 2
        start = midpoint - args.subtract_from_start
        end = midpoint + args.add_to_end
    if start < 1:
        out.close()
        stop_err('Requested expansion places region beyond chromosome bounds.')
    new_line = '\t'.join([items[0], items[1], items[2], str(start), str(end), items[5], items[6], items[7], items[8]])
    out.write(new_line)
out.close()

