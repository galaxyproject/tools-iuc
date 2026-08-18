#!/usr/bin/python

import sys
import math
import argparse

def main():
    parser = argparse.ArgumentParser(prog='blastfiler', usage='%(prog)s [options]', description='Parses blast output for coverage and pident [static (pident and length) or rost1999] to get subset of input. The options "covs" and "covq" describes the coverage per subject and per query by dividing the alignment length by the length of the subject or the query. Here one has to choose the correct "type" option so that the length is assigned to aa or bp. E.g. if  blastp was used, both qlen and slen  will be in aa. If blastx was used, qlen will be in bp and slen in aa. If tblastn was used qlen will be in aa and slen will be in bp. If tbalstx was used qlen and slen will be in bp. If blastn was used qlen and slen will be in bp. Depending on the blast command which was used to generate the output one can choose bwetween the following options for the "outfmt" option:\nstd: -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"\nraw: -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen nident gaps score"\nrbhplus:  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen nident gaps score"')
    parser.add_argument('-i', help='specify blast input file')
    parser.add_argument('-o', help='specify blast output file')
    parser.add_argument('-p', default=float(30), type=float, help='specify pident cutoff [default: 30]')
    parser.add_argument('-l', default=float(80), type=float, help='specify alignment length cutoff [default: 80]')
    parser.add_argument('-pid', default='static', choices=['static', 'rost1999'], help='choose pident filtering either as static [pident+length] or rost1999 [default: static]')
    parser.add_argument('-covs', default=float(0), type=float, help='specify covs threshol (length/slen), [default: 0.0]')
    parser.add_argument('-covq', default=float(0), type=float, help='specify covq threshol (length/qlen), [default: 0.0]')
    parser.add_argument('-t', default='blastp', help='specify type of blast which generated the output table, [default: "blastp"]')
    parser.add_argument('-outfmt', default='std', help='specify outfmt which was used with the blast command to generate the output table, [default: "std"]')
    args = parser.parse_args()
    if args.i is None:
        parser.print_help()
        sys.exit('\nPlease specify input file')
    species_A_infile = args.i
    if args.o is None:
        parser.print_help()
        sys.exit('\nPlease specify output file')
    outfile = args.o
    option_pident = args.p
    option_length = args.l
    pident_method = args.pid
    covs_option = args.covs
    covq_option = args.covq
    outfmt = args.outfmt
    t = args.t
    
    # print("command arguments used:")
    # print(args)
    print("Start blastpfilter")
    
    def get_pident_by_length(x):
        if x<=11:
            return float(100)
        if x<=450:
            return 480*float(x)**(-0.32*(1+math.exp(-float(x)/float(1000))))
        if x>450:
            return 19.5
    
    ##species A
    blast_species_A_dict={}
    blast_species_A = open(species_A_infile,'r')
    
    with open(outfile,'w') as handle:
        for line in blast_species_A:
            parts = line.strip().split('\t')
            if outfmt=='raw':
                qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qcovs,qcovhsp,qlen,slen,nident,gaps,score = parts
            if outfmt=='rbhplus':
                qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen,nident,gaps,score = parts
            if outfmt=='std':
                qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore = parts
            evalue = float(evalue)
            bitscore = float(bitscore)
            pident = float(pident)
            length = float(length)
            if outfmt=='raw':
                qcovs = float(qcovs)
                qcovhsp = float(qcovhsp)
            if outfmt=='raw' or outfmt=='rbhplus':
                if t=='blastp' or t=='blastn':
                    qlen = float(qlen)
                    slen = float(slen)
                if t=='blastx':
                    qlen = float(qlen)/3
                    slen = float(slen)
                if t=='tblastn':
                    qlen = float(qlen)
                    slen = float(slen)/3
                if t=='tblastx':
                    qlen = float(qlen)/3
                    slen = float(slen)/3
                nident = float(nident)
                gaps = float(gaps)
                score = float(score)
                covs = length/slen
                covS = qlen/slen
                covns = nident/slen
                covq = length/qlen
                covQ = slen/qlen
                covnq = nident/qlen
            if pident_method=='static':
                if outfmt=='raw' or outfmt=='rbhplus':
                    if pident>=option_pident and length>=option_length and covs>=covs_option and covq>=covq_option:
                        handle.write(line)
                if outfmt=='std':
                    if pident>=option_pident and length>=option_length:
                        handle.write(line)
            if pident_method=='rost1999':
                if outfmt=='raw' or outfmt=='rbhplus':
                    if pident>=get_pident_by_length(length) and covs>=covs_option and covq>=covq_option:
                        handle.write(line)
                if outfmt=='std':
                    if pident>=get_pident_by_length(length):
                        handle.write(line)
    print("Finished blastpfilter")

if __name__ == '__main__':
    main()
