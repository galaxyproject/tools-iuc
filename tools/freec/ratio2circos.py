import argparse
import os
import math

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, default='./output/sample.bam_ratio.BedGraph', type=str)
parser.add_argument('-o', '--output', required=True, default='./output/sample.bam_ratio_log2_circos.txt', type=str)
parser.add_argument('-p', '--ploidy', required=True, default=2, type=int)
args = parser.parse_args()

path = os.path.dirname(args.input)
output = os.path.join(path, args.output)

with open(args.input) as file:
    for line in file.readlines():
        ls = line.split()
        if ls[0] != "track" and float(ls[3]) > 0:
            log2_ratio = math.log2(float(ls[3]) / args.ploidy)
            with open(output, "a") as out:
                out.write("{}\t{}\t{}\t{}\n".format(ls[0], ls[1], ls[2], log2_ratio))
