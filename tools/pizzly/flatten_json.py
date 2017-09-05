import json
import sys
from collections import OrderedDict


def loadJSON(fn):
    with open(fn) as f:
        JJ = json.load(f, object_pairs_hook=OrderedDict)
    return JJ['genes']


def outputGeneTable(fusions, outf, filters=None):
    outf.write('\t'.join("geneA.name geneA.id geneB.name geneB.id paircount splitcount transcripts.list".split()))
    outf.write('\n')
    for gf in fusions:
        gAname = gf['geneA']['name']
        gAid = gf['geneA']['id']
        gBname = gf['geneB']['name']
        gBid = gf['geneB']['id']
        pairs = str(gf['paircount'])
        split = str(gf['splitcount'])
        txp = [tp['fasta_record'] for tp in gf['transcripts']]

        outf.write('\t'.join([gAname, gAid, gBname, gBid, pairs, split, ';'.join(txp)]))
        outf.write('\n')


def usage():
    print("Usage: python flatten_json.py fusion.out.json [genetable.txt]")
    print("")
    print("       outputs a flat table listing all gene fusions, if the output file is not")
    print("       specified it prints to standard output")


if __name__ == "__main__":
    nargs = len(sys.argv)
    if nargs <= 1:
        usage()
    else:
        infn = sys.argv[1]
        fusions = loadJSON(infn)
        outf = sys.stdout
        if nargs == 3:
            outf = open(sys.argv[2], 'w')

        outputGeneTable(fusions, outf)

        if outf != sys.stdout:
            outf.close()
