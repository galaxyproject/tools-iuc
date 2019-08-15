#!/usr/bin/env python
import csv
import os
import sys


def get_description(line, parent):
    cols = line.split('\t')
    labels = {}
    for pair in cols[8].split(";"):
        k, v = pair.split('=')
        labels[k] = v

    if (cols[2]) == "CDS" and labels["Parent"] == parent:
        return labels.get("Note", '-')
    return '-'


def convert_to_prot_table(fname, output_name):
    gff_file = open(fname)
    output_file = open(output_name, 'w')
    writer = csv.writer(output_file, delimiter='\t')
    lines = gff_file.readlines()
    gff_file.close()
    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith('#'):
            continue
        cols = line.split('\t')
        if (len(cols) < 9):
            print("Ignoring invalid row with entries: {0}".format(cols))
        elif (cols[2]) == "region":
            continue
        elif (cols[2]) == "CDS":
            continue
        elif (cols[2]) == "gene":
            start = int(cols[3])
            end = int(cols[4])
            strand = cols[6].strip()
            labels = {}
            diff = int(abs(end - start) / 3)  # What is this called?
            for pair in cols[8].split(";"):
                k, v = pair.split('=')
                labels[k.strip()] = v.strip()

            Rv = labels["locus_tag"].strip()  # error out if not found
            gene = labels.get('Name', '')
            desc = get_description(lines[i + 1], labels.get("ID", "")) if (i + 1) < len(lines) else '-'
            vals = [desc, start, end, strand, diff, '-', '-', gene, Rv, '-']
            writer.writerow(vals)
    output_file.close()


if __name__ == "__main__":
    usage_string = "Usage: python gff-prot-converter.py <gff filename> <output filename>"

    if len(sys.argv) < 3:
        print(usage_string)
        sys.exit(0)
    file_name = sys.argv[1]
    if not os.path.exists(file_name):
        print("File not found. Exiting...")
        print(usage_string)
        sys.exit(0)
    convert_to_prot_table(file_name, sys.argv[2])
