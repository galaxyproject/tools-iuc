#!/usr/bin/env python

# filter names and nodes dmp files by a list of given IDs
# parent node IDs will be added if needed
#
# IDs will be renamed to give a consecuive set of IDs: 1,2,...
# oderwise dmnd databases including taxonomy will be huge
# also make make sure that the order of the taxids is not changed

from sys import argv

names_file_name = argv[1]
nodes_file_name = argv[2]
prot2ids_file_name = argv[3]
names_file_out_name = argv[4]
nodes_file_out_name = argv[5]
prot2ids_file_out_name = argv[6]

parent = dict()
with open(nodes_file_name) as nodes_file:
    for line in nodes_file:
        line = line.strip().split("|")
        parent[line[0].strip()] = line[1].strip()

initial_ids = set()
with open(prot2ids_file_name) as prot2ids_file:
    for i, line in enumerate(prot2ids_file):
        if i == 0:
            continue
        line = line.strip().split()
        initial_ids.add(line[2].strip())

ids = set()
while len(initial_ids):
    i = initial_ids.pop()
    p = parent[i]
    if p == i:
        ids.add(p)
        continue
    ids.add(i)
    initial_ids.add(p)

id_map = dict()
with open(names_file_name) as names_file, open(names_file_out_name, "w") as names_file_out:
    for line in names_file:
        line = line.strip().split("|")
        id = line[0].strip()
        if id not in ids:
            continue
        if id not in id_map:
            id_map[id] = len(id_map) + 1
        names_file_out.write(f'{id_map[id]}\t|{"|".join(line[1:])}\n')

print(f'taxonlist for test 2 needs to be {id_map["33090"]}')

with open(nodes_file_name) as nodes_file, open(nodes_file_out_name, "w") as nodes_file_out:
    for line in nodes_file:
        line = line.strip().split("|")
        node = line[0].strip()
        parent = line[1].strip()
        if node not in ids or parent not in ids:
            continue
        nodes_file_out.write(f'{id_map[node]}\t|\t{id_map[parent]}\t|{"|".join(line[2:])}\n')

with open(prot2ids_file_name) as prot2ids_file, open(prot2ids_file_out_name, "w") as prot2ids_file_out:
    for i, line in enumerate(prot2ids_file):
        if i == 0:
            prot2ids_file_out.write(line)
            continue
        line = line.strip().split()
        id = line[2].strip()
        line[2] = str(id_map[id])
        prot2ids_file_out.write("\t".join(line) + "\n")
