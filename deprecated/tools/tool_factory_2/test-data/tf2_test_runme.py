# bog standard argparse for 3 possible comma separated lists
# followed by some silly reverse each row code provided as an example
# you're supposed to replace it with your great code..
import argparse
import copy

argp = argparse.ArgumentParser()
argp.add_argument('--INNAMES', default=None)
argp.add_argument('--INPATHS', default=None)
argp.add_argument('--OUTPATH', default=None)
argp.add_argument('--additional_parameters', default=[], action="append")
argp.add_argument('otherargs', nargs=argparse.REMAINDER)
args = argp.parse_args()
fout = open(args.OUTPATH, 'w')
sins = open(args.INPATHS.split(',')[0]).readlines()
for i, sin in enumerate(sins):
    row = sin.strip().split('\t')
    rrow = copy.copy(row)
    lrow = len(row)
    if (lrow > 1):
        for j in range(lrow):
            rrow[j] = row[lrow - j - 1]
        fout.write('\t'.join(rrow))
        fout.write('\n')
fout.close()
