import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--folder', default='chromosomes', type=str)
parser.add_argument('-p', '--prefix', default='chr', type=str)
parser.add_argument('-s', '--suffix', default='.fa', type=str)
args = parser.parse_args()

f_content = os.listdir(args.folder)

if f_content[0][:1] == ">":
    for file in f_content:
        chr_nr = file.split()[0][1:]
        os.rename(args.folder+"/"+file, args.folder+"/"+args.prefix+chr_nr+args.suffix)
