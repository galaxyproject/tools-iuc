# Try to convert constellations files into the format expected by lineagespot.
# Constellations files can define parent lineages, in which case the script
# parses parent mutations recursively and adds them to the signature of the
# child.

# CURRENT AND GENERAL LIMITATIONS
# Important to understand, please read carefully
# 1. Constellations sometimes uses base instead of amino acid positions for
# defining mutations. These can take two forms like in these examples:
# "nuc:C8986T", i.e. a SNV given in base coordinates
# "del:22029:6", i.e. a deletion of 6 bases given in base coordinates
# The current version of the script makes no attempt to convert such lines to
# amino acid coordinates, but simply drops them.
# 2. In other cases, constellations lists deletions in amino acid poisitions like
# this:
# "s:HV69-"
# While this notation could be parsed such lines are currently *also* dropped
# because it's not entirely clear how lineagespot describes deletions.
# 3. In some cases, constellation also provides mutations in mature peptide
# coordinates, like "nsp15:K259R". Lines like this are currently dropped, too.
# 4. The constellations data provided by
# https://github.com/cov-lineages/constellations
# lists mostly lineage-defining mutations that can be used to *distinguish*
# between lineages, but makes no attempt to provide complete lists of mutations
# (even through parent lineage definitions) for any lineage.

import argparse
import json
import os
import re
import sys


genes_names_translation = {
    "orf1a": "ORF1a",
    "orf1ab": "ORF1ab",
    "1ab": "ORF1ab",
    "orf1b": "ORF1b",
    "s": "S",
    "spike": "S",
    "orf3a": "ORF3a",
    "e": "E",
    "m": "M",
    "n": "N",
    "orf6": "ORF6",
    "orf7a": "ORF7a",
    "orf7b": "ORF7b",
    "orf8": "ORF8",
    "8": "ORF8",
    "n": "N",
    # NOTE: in constellations, mutations are sometimes, but not always, given
    # in nsp coordinates instead of ORF1a/b ones. Currently, we drop these,
    # while we should convert instead!!!
    "nsp2": "NSP2",
    "nsp3": "NSP3",
    "nsp4": "NSP4",
    "nsp5": "NSP5",
    "nsp6": "NSP6",
    "nsp7": "NSP7",
    "nsp8": "NSP8",
    "nsp9": "NSP9",
    "nsp10": "NSP10",
    "nsp12": "NSP12",
    "nsp13": "NSP13",
    "nsp14": "NSP14",
    "nsp15": "NSP15",
    "nsp16": "NSP16",
}


lineagespot_template = dict.fromkeys(["ORF1a", "ORF1b", "S", "ORF3a", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "E", "N"])
definitions = {}

pat = re.compile(r'(?P<gene>.+):(?P<ref>[A-Z]+)(?P<pos>\d+)(?P<alt>[A-Z*]+)')


def read_lineage_variants(x, lineage_name):
    data = json.load(x)

    sites = {}
    for mut in data["sites"]:
        match = pat.match(mut)
        if match is None:
            # Likely a del or nuc mutation given at the base level
            continue
        # try to get a canonical gene name
        gene = genes_names_translation.get(
            match.group('gene'),
            match.group('gene')
        )
        pos = int(match.group('pos'))
        if gene == 'ORF1ab':
            # constellations isn't very consistent in representing ORF1ab
            # mutations. They may be provided in ORF1a or ORF1b coordinates,
            # but could also just be given as ORF1ab.
            if pos <= 4401:
                gene = 'ORF1a'
            else:
                gene = 'ORF1b'
                # 4715 == 314 in constellations
                pos = pos - 4401
        if gene not in sites:
            sites[gene] = {}
        sites[gene][pos] = (match.group('ref'), match.group('alt'))

    # recursively parse parent lineages and
    # add their mutations to the global definitions
    if "parent_lineage" in data["variant"]:
        x_parent = data["variant"]["parent_lineage"]
        if x_parent not in definitions:
            parent_filename = f"c{x_parent}.json"
            lineage_def_dir = os.path.dirname(x.name)
            parent_file = os.path.join(lineage_def_dir, parent_filename)
            if not os.path.isfile(parent_file):
                raise FileNotFoundError(
                    f"{x_parent} is defined as a parent of {lineage_name}, but "
                    f"definitions file {parent_filename} not found in "
                    f"{lineage_def_dir}."
                )
            with open(parent_file) as parent_in:
                read_lineage_variants(parent_in, x_parent)

        # update the sites dictionary to include also mutations defined for the parent
        for gene, muts in definitions[x_parent].items():
            if gene in sites:
                for pos, ref_alt in muts.items():
                    if pos in sites[gene]:
                        # exotic case of a parent site being affected in the child
                        # lineage again. Kepp the child site unaltered.
                        continue
                    sites[gene][pos] = ref_alt
            else:
                # only the parent has mutations in this gene listed
                sites[gene] = muts
    # done with this lineage and all of its parents
    definitions[lineage_name] = sites


parser = argparse.ArgumentParser()
parser.add_argument(
    "-i", "--input", required=True,
    help="Name of the input folder"
)
parser.add_argument(
    "-o", "--output", required=True,
    help="Name of the output folder"
)
if len(sys.argv) < 2:
    sys.exit('Please run with -h / --help for help.')

args = parser.parse_args()

for definitions_file in os.listdir(args.input):
    # In constellations, the only reliable way to get the lineage name is from
    # the file name by stripping the .json suffix from it and dropping the
    # leading 'c' (e.g. cBA.5.json holds the definition for lineage BA.5).
    if definitions_file[0] != 'c' or definitions_file[-5:] != '.json':
        continue
    lineage_name_from_file = definitions_file[1:-5]
    if lineage_name_from_file in definitions:
        # seems we have parsed this lineage already as a parent of another lineage
        continue
    with open(os.path.join(args.input, definitions_file)) as data_in:
        read_lineage_variants(data_in, lineage_name_from_file)

for lineage, sites in definitions.items():
    # if path isn't there, create one could be added
    with open(os.path.join(args.output, lineage) + '.txt', "w") as data_out:
        data_out.write('gene\tamino acid\n')
        for gene, muts in sites.items():
            if gene in lineagespot_template:
                for pos, ref_alt in muts.items():
                    data_out.write(f'{gene}\t{ref_alt[0]}{pos}{ref_alt[1]}\n')
