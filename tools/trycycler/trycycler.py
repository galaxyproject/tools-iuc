#!/usr/bin/python3

from os import walk,listdir,path, makedirs,remove
from sys import argv
from shutil import copy2
import tarfile


def cluster(output_folder):
    counter = 1
    output_folder = argv[2]
    for root, dir, files in walk(output_folder):
        if  "1_contigs" in dir:
            full_path = path.join(root,dir[0])
            files = [path.join(full_path,x) for x in listdir(full_path)]
            with open(path.join(output_folder,"cluster_0{}.fasta".format(counter)),"a") as tmp:
                [tmp.write(open(x).read()) for x in files]
            counter+=1
            continue


def reconcile(input_file):
    number_cluster = [x for x in input_file[-1:0:-1] if x.isdigit()][0]
    fullpath="selected_cluster/cluster_0{}/1_contigs/".format(number_cluster)    
    with open(input_file) as tmp:
        reads = [">"+x for x in (tmp.read().split(">"))[1:]]
        for read in reads:
            fasta_name = read.split("\n")[0][1:]+".fasta"
            output_fasta = "{}{}".format(fullpath,fasta_name)
            with open(output_fasta,"w") as temporal:
                temporal.write(read)
            
def main():
    if argv[1] == "cluster": cluster(argv[2])
    if argv[1] == "reconcile": reconcile(argv[2])

if __name__ == "__main__":
    main()
    
