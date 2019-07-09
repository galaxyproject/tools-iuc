#!/usr/bin/env python3


import sys
import argparse


def parse_args():
    parser = argparse.ArgumentParser(prog='deneme.py')
    hg19_flag = False
    uninform = False
    zeromaf = False
    usehapmap = True
    interval = max_interval = centim = -1

    parser.add_argument('mapFile', nargs=1, help='Map File')
    parser.add_argument('mafFile', nargs=1, help='Maf File')
    parser.add_argument('genotypesFile', nargs=1, help='Genotypes File')
    parser.add_argument('--hg19', default=True, action='store_true', help='use human genome ver 19 (instead of the version 18)')
    parser.add_argument('--nohapmap', default=False, action='store_true', help='to not check hapmap data')
    parser.add_argument('-u', '--uninformative', default=True, action='store_true', help='to include uninformative markers too')
    parser.add_argument('-z', '--zeromaf', default=True, action='store_true', help='to include zero minor allele frequency markers')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-s', '--spacing', type=int, help='select markers at intervals of N, or specify the max number of markers N, or select at intervals of N centiMorgans. Both -u and -z are autoselected for this flag.')
    group.add_argument('-m', '--max', type=int, help='select markers at intervals of N, or specify the max number of markers N, or select at intervals of N centiMorgans. Both -u and -z are autoselected for this flag.')
    group.add_argument('-cm', '--centim', '--centimorgans', type=float, help='select markers at intervals of N, or specify the max number of markers N, or select at intervals of N centiMorgans. Both -u and -z are autoselected for this flag.')

    args = vars(parser.parse_args())
    interval = args['spacing']
    usehapmap = args['nohapmap']
    hg19_flag = args['hg19']
    uninform = args['uninformative']
    zeromaf = args['zeromaf']
    max_interval = args['max']
    centim = args['centim']
    files = [args['mapFile'][0], args['mafFile'][0], args['genotypesFile'][0]]
    return files, interval, max_interval, hg19_flag, uninform, zeromaf, centim, usehapmap


# =====================================================


def uninformativeSNPs(genotypesfile):
    uninformative = {}
    genotypecount = 0

    print("[INFO] reading '%s'..." % genotypesfile, file=sys.stderr)

    f = open(genotypesfile)
    for line in f:
        # if line.startswith('rs'):
        tmp = line.strip().split()
        marker = tmp[0]
        gtypes = tmp[1:]
        gtypes = set([x for x in gtypes if x != 'NC'])
        uninformative[marker] = True if len(gtypes) == 1 else False
        # genotypecount += 1
        genotypecount += 1
    f.close()

    print("[INFO] found %d/%d SNPs uninformative" % (len([x for x in list(uninformative.values()) if x is True]), genotypecount), file=sys.stderr)
    return uninformative

# ==================


def minorAlleleFreqs(maffile):
    frequencies = {}

    print("[INFO] reading '%s'..." % maffile, file=sys.stderr)

    f = open(maffile)
    f.readline()  # strip header
    for line in f:
        # if line.startswith('rs'):
        if len(line) < 5:
            continue
        # try:
        rs, freq = line.strip().split()
        frequencies[rs] = float(freq)
    # except ValueError:
    # print >> sys.stderr, '\n', line
    # exit(-1)
    f.close()

    print("[INFO] read %d SNP minor allele frequencies" % len(frequencies), file=sys.stderr)

    return frequencies

# ==================


def get_hapmap_map(i):
    data = {}
    # print >> sys.stderr, "[INFO] reading c%d reference..." % i
    global hg19_flag

    try:
        i = int(i)
    except ValueError:
        pass

    i = str(i)

    try:
        if hg19_flag:
            g = open(hapmap_path + "genetic_map_GRCh37_chr%s.txt" % i)
        else:
            g = open(hapmap_path + "genetic_map_chr%s_b36.txt" % i)
    except IOError:
        # error_msg = "No Hapmap data for chr"+i+(" hg19" if hg19_flag else " hg18")+" - using uncorrected cM positions in original map"
        # print >> sys.stderr, "\n[Caution]",error_msg,
        # print >> logfile, "error: "+error_msg    #log message in readme_inp

        return -1

    g.readline()

    if hg19_flag:
        for line in g:
            chr_nm, hm_phy, hm_com, hm_gen = line.strip().split()
            data[hm_phy] = float(hm_gen)

    else:
        for line in g:
            hm_phy, hm_com, hm_gen = line.strip().split()
            data[hm_phy] = float(hm_gen)

    g.close()
    return data

# ==================

# Default method of finding all the informative SNPs in the map file


def informativeSNPs(mapfile, uninformative, frequencies, spacing, zeromaf, use_uninf, checkhap, interval=0.05):
    geneticmap_data = None
    markercount = 0
    threshold = 0.0

    print("[INFO] reading '%s' and printing new map..." % mapfile, file=sys.stderr)
    print("[INFO] chromosome", end=' ', file=sys.stderr)

    f = open(mapfile)

    print(f.readline().strip())
    current_chr = 0
    line_count = 0

    for line in f:
        chr, rsid, genpos, physpos, rsid2, junk = line.strip().split()
        genpos = float(genpos)
        # try:
        # chr = int(chr)
        # except:
        # pass # X, XY, Y
        # print line.strip()
        # continue

        if chr != current_chr:
            current_chr = chr
            threshold = 0.0
            print("%s" % chr, end=' ', file=sys.stderr)
            geneticmap_data = get_hapmap_map(chr)

        skip = False

        if checkhap:
            if geneticmap_data != -1:  # Has hapmap data:
                if physpos in geneticmap_data:
                    genpos = geneticmap_data[physpos]  # Update genpos with hapmap data
                else:
                    skip = True  # Has hapmap, no SNP? Skip

        if not(skip):
            if spacing != -1:
                if line_count % spacing == 0:
                    print('\t'.join([chr, rsid, str(genpos), physpos, rsid, junk]))
                    markercount += 1
            else:
                if genpos > threshold:
                    try:
                        # if not rsid.startswith("rs") :  # not kgp friendly this
                        # continue
                        if not(use_uninf):
                            if uninformative[rsid]:
                                continue

                        if not(zeromaf):
                            if (frequencies[rsid] == 0) or (frequencies[rsid] == 1):
                                continue

                    except KeyError:
                        continue

                    threshold = genpos + interval

                    print('\t'.join([chr, rsid, str(genpos), physpos, rsid, junk]))
                    markercount += 1

        line_count += 1

    f.close()
    print("\n[INFO] new map contains %d SNPs" % markercount, file=sys.stderr)


# MAIN

if __name__ == "__main__":
    files, equally_spaced, max_interval, hg19_flag, use_uninf, zeromaf, cm, checkhap = parse_args()

    mapfile = files[0]
    maffile = files[1]
    genotypesfile = files[2]

    hapmap_path = '/usr/local/maps/hapmap_geneticmaps/'

    frequencies = [] if zeromaf else minorAlleleFreqs(maffile)

    uninformative = [] if use_uninf else uninformativeSNPs(genotypesfile)
    # print >> sys.stderr, uninformative

    # Parse Spacing args
    spacing = -1

    if equally_spaced is not None:
        spacing = equally_spaced
        spacing = 1 if (spacing < 1) else spacing

    if max_interval is not None:
        def num_lines_mapf(fname):
            i = 0
            with open(fname) as f:
                for i, l in enumerate(f):
                    pass
            return i + 1
        spacing = int( num_lines_mapf(mapfile) / float(max_interval) ) - 1
        spacing = 1 if (spacing < 1) else spacing

    # Parse centimorgan spacing
    interval = 0.05 if (cm is None) else cm
    mark = informativeSNPs(mapfile, uninformative, frequencies, spacing, zeromaf, use_uninf, checkhap, interval )
