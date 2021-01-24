#!/usr/bin/env python3

from os import path, walk
from sys import argv


def cluster(output_folder):
    counter = 1
    for root, dir, files in walk(output_folder):
        if root.endswith('1_contigs'):
            output_path = path.join(output_folder , f"cluster_0{counter}.fasta")
            with open(output_path, "a") as out_cluster:
                for fasta in files:
                    fasta_path = path.join(root, fasta)
                    fasta = open(fasta_path).read()
                    out_cluster.write(fasta)
            counter += 1


def reconcile(input_file):
    number_cluster = [x for x in input_file[-1:0:-1] if x.isdigit()][0]
    full_path = f"selected_cluster/cluster_0{number_cluster}/1_contigs/"    
    with open(input_file) as tmp:
        for line in tmp:
            if ">" in line:
                filename = line[1:].strip()
                output_fasta = f"{full_path}{filename}.fasta"
            open(output_fasta , "a").write(line)


def main():
    if argv[1] == "cluster": cluster(argv[2])
    if argv[1] == "reconcile": reconcile(argv[2])


if __name__ == "__main__":
    main()
