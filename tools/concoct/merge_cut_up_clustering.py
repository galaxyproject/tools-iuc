#!/usr/bin/env python

import argparse
import re
import sys
from collections import Counter
from collections import defaultdict


CONTIG_PART_EXPR = re.compile(r'(.*)\.concoct_part_([0-9]*)')


def original_contig_name_special(contig_id):
    try:
        original_id, part_index = CONTIG_PART_EXPR.match(contig_id).group(1, 2)
        return original_id, part_index
    except AttributeError:
        # No matches for concoct_part regex.
        return contig_id, 0


parser = argparse.ArgumentParser()
parser.add_argument("--input", action="store", dest="input", help="Tabular file with cut up clusters")
parser.add_argument("--output", action="store", dest="output", help="Output file with merged clusters")

args = parser.parse_args()

# Get cut up clusters
all_seqs = {}
all_originals = defaultdict(dict)
with open(args.input, 'r') as ifh:
    for i, line in enumerate(ifh):
        if i == 0:
            if 'contig_id' not in line:
                sys.stderr.write("ERROR nvalid clustering file, 'contig_id' is not found in the header.")
                sys.exit(-1)
            # Skip header.
            continue
        line = line.rstrip('\r\n')
        contig_id, cluster_id = line.split('\t')
        original_contig_name, part_id = original_contig_name_special(contig_id)
        all_originals[original_contig_name][part_id] = cluster_id

# Merge cut up clusters.
with open(args.output, 'w') as ofh:
    ofh.write("contig_id\tcluster_id\n")
    for original_contig_id, part_ids_d in all_originals.items():
        if len(part_ids_d) > 1:
            c = Counter(part_ids_d.values())
            cluster_id = c.most_common(1)[0][0]
            c_string = [(a, b) for a, b in c.items()]
            # Here if len(c.values()) > 1,
            # then no cluster for contig.
        else:
            cluster_id = list(part_ids_d.values())[0]
        ofh.write(f"{original_contig_id}\t{cluster_id}\n")
