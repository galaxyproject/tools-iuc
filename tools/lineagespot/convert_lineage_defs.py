import json
import os
import sys
import argparse

genes_names_translation = {
               "orf1a": "ORF1a",
               "ORF1a": "ORF1a",
               "orf1ab": "ORF1ab",
               "ORF1ab": "ORF1ab",
               "1ab": "ORF1ab",
               "orf1b": "ORF1b",
               "ORF1b": "ORF1b",
               "s": "S",
               "S": "S",
               "spike": "S" ,
               "orf3a": "ORF3a",
               "ORF3a": "ORF3a",
               "M": "M",
               "m": "M",
               "N": "N",
               "n": "N",
               "orf6": "ORF6",
               "ORF8": "ORF8",
               "orf8": "ORF8",
               "8": "ORF8",
               "ORF7a": "ORF7a",
               "ORF7b": "ORF7b",
               "XE-parent1": "XE-parent1",
               "E": "E",
               "e": "E",
               "NSP2": "NSP2",
               "nsp3": "NSP3",
               "nsp4": "NSP4",
               "nsp5": "NSP5",
               "nsp6": "NSP6",
               "nsp7": "NSP7",
               "nsp12": "NSP12",
               "nsp13": "NSP13",
               "nsp15": "NSP15"}
lineagespot_template = dict.fromkeys(["ORF1a", "ORF1b", "S", "ORF3a", "M", "ORF7a", "ORF7b", "ORF8", "E", "N"])
definitions = {}
def read_lineage_variants(x):

# opening
    data = json.load(x)

    label_initial = lineagespot_file
    if label_initial not in definitions:

# initial appending
        sites = {}
        for genes in data["sites"]:
            if genes[:3] == "del" or genes[:3] == "nuc":
                    continue
            k1, v1 = genes.split(":")
            translated_gene = genes_names_translation.get(k1)             # comment
            if translated_gene not in sites:
                sites[translated_gene] = []
            sites[translated_gene].append(v1)
        definitions[label_initial] = sites

# parent opening
        if "parent_lineage" in data["variant"]:
            x_parent = data["variant"]["parent_lineage"]
            if x_parent not in definitions:
                f_parent = open("/home/eliberg/constellations/constellations/definitions/c"+x_parent+".json")
                read_lineage_variants(f_parent)

    return



#parser
parser = argparse.ArgumentParser()

#parser.add_argument("-h" "--help")

parser.add_argument(
    "-i", "--input", required=True,
    help="Name of the Input folder"
)
parser.add_argument(
    "-o", "--output", required=True,
    help="Name of the output folder"
)
if len(sys.argv) < 2:
    print(notice, file=sys.stderr)
    print('\nPlease run with -h / --help for help.', file=sys.stderr)
    sys.exit(2)

args = parser.parse_args()

input = args.input
for definitions_file in os.listdir(input):
    lineagespot_file = definitions_file[:-5]
    definitions_file_path = os.path.join(args.input,  definitions_file)
    open_definitions = open(definitions_file_path)
    read_lineage_variants(open_definitions)

output = args.output
for lineage, genes in definitions.items():
    # if path isn't there, create one could be added
    with open(os.path.join(output, lineage), "w") as file:
        file.write("%s:%s\n" % ("gene", "amino acid"))
        for gene in genes:
            if gene in lineagespot_template:
                for gene_site in genes[gene]:
                    file.write("%s:%s\n" % (gene, gene_site))

