import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genome', type=str)
parser.add_argument('-f', '--folder', default='chromosomes', type=str)
parser.add_argument('-p', '--prefix', default='chr', type=str)
parser.add_argument('-s', '--suffix', default='.fa', type=str)
args = parser.parse_args()

os.makedirs(args.folder, exist_ok=True)

if args.genome:
    with open(args.genome, "r") as fasta:
        for line in fasta:
            if line[:1] == ">":
                chr_nr = line[1:].split()[0]
            with open(args.folder + "/" + args.prefix + chr_nr + args.suffix, "a") as chr_fasta:
                chr_fasta.write(line)
