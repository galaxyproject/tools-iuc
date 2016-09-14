#!/usr/bin/env python
import argparse
import subprocess
import os
import shlex

def createIndex(fname):
    """
    Create an index for a file, return where it is
    """
    d = os.mkdir("index_dir")
    os.link(fname, "index_dir/genome.fa")
    cmd = "{}/bwameth.py index {}/index_dir/genome.fa".format(os.path.dirname(__file__), os.getcwd())
    proc = subprocess.Popen(args=shlex.split(cmd))
    returncode = proc.wait()
    if returncode != 0:
        raise Exception("Error during '%s'" % cmd)
    return "{}/index_dir/genome.fa".format(os.getcwd())


parser = argparse.ArgumentParser(description="A wrapper around bwameth for Galaxy. There's minimal argument checking done")
parser.add_argument('-p', '--numThreads', type=int, default=4, help="number of threads")
parser.add_argument('--makeIndex', help="Given a fasta file, index it")
parser.add_argument('--premadeIndex', help="If an index already exists, this is the fasta file associated with it (the index files must be in the same directory")
parser.add_argument('--readGroup', help="The read group text, if desired")
parser.add_argument('files', nargs="+", help="Fasta files (possibly gzipped, if the file names end in .gz)")
args = parser.parse_args()

if args.makeIndex:
    ifile = createIndex(args.makeIndex)
else:
    ifile = args.premadeIndex

files = " ".join(['{}'.format(x) for x in args.files])
if args.readGroup:
    files = "{} --read-group '{}'".format(files, args.readGroup)

cmd = "{}/bwameth.py -t {} --reference '{}' {} | samtools view -u - | samtools sort -@ {} - output".format(os.path.dirname(__file__), args.numThreads, ifile, files, args.numThreads)
of = open("temp.sh", "w")
of.write("{}\n".format(cmd))
of.close()
try:
    subprocess.check_call(["bash", "./temp.sh"])
except:
    raise Exception("Error during '%s'" % cmd)
