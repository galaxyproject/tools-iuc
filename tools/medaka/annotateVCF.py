#!/usr/bin/env python3

### Takes in VCF file and a samtools mpileup output file
### Fills in annotation for the VCF file including AF, DP
### SB, and DP4
###
### Usage statement:
### python annotateVCF.py in_vcf.vcf in_mpileup.txt out_vcf.vcf
###
### Can generate in_mileup.txt with samtools mpileup (and can restrict which sites to generate pileups for with in_vcf.vcf)

### 06/18/2020 - Nathan P. Roach
### Written in the employ of GalaxyWorks LLC

import sys
from math import log10,isnan
from scipy.stats import fisher_exact

def pval_to_phredqual(pval):
    return int(round(-10. * log10(pval)))

def parseSimpleSNPpileup(fields,ref_base,alt_base):
    base_to_idx = {
        'A' : 0,
        'a' : 0,
        'T' : 1,
        't' : 1,
        'C' : 2,
        'c' : 2,
        'G' : 3,
        'g' : 3
    }

    base_to_idx_stranded = {
        'A' : 0,
        'T' : 1,
        'C' : 2,
        'G' : 3,
        'a' : 4,
        't' : 5,
        'c' : 6,
        'g' : 7
    }
    ref_base2 = fields[2]
    counts = [0,0,0,0]
    stranded_counts = [0,0,0,0,0,0,0,0]
    ref_idx = base_to_idx[fields[2]]
    dp = int(fields[3])
    carrot_flag = False
    ins_flag = False
    ins_str = ""
    ins_len = 0
    insertion = ""
    del_flag = False
    del_str = ""
    del_len = 0
    deletion = ""
    # dollar_flag = False
    for base in fields[4]:
        if carrot_flag:
            carrot_flag = False
            continue
        if ins_len > 0:
            insertion += base
            ins_len -= 1
            continue
        if del_len > 0:
            deletion += base
            del_len -= 1
            continue
        if ins_flag:
            if base.isdigit():
                ins_str += base
            else:
                ins_len = int(ins_str) - 1
                insertion = base
                ins_flag = False
        elif del_flag:
            if base.isdigit():
                del_str += base
            else:
                del_len = int(del_str) - 1
                deletion = base
                del_flag = False
        else:
            if base == '^':
                carrot_flag = True
                continue
            elif base == '$':
                continue
            elif base == '+':
                ins_flag = True
            elif base == '-':
                del_flag = True
            elif base == '.':
                counts[ref_idx] += 1
                stranded_counts[base_to_idx_stranded[ref_base2]] += 1
            elif base == ',':
                counts[ref_idx] += 1
                stranded_counts[base_to_idx_stranded[ref_base2.lower()]] += 1
            elif base == 'N' or base == 'n':
                continue
            elif base == '*':
                continue
            else:
                counts[base_to_idx[base]] += 1
                stranded_counts[base_to_idx_stranded[base]] += 1
    af = float(counts[base_to_idx[alt_base]]) / float(sum(counts))
    if float(sum(stranded_counts[0:4])) == 0:
        faf = float("nan")
    else:
        faf = float(stranded_counts[base_to_idx_stranded[alt_base]]) / float(sum(stranded_counts[0:4]))
    if float(sum(stranded_counts[4:])) == 0:
        raf = float("nan")
    else:
        raf = float(stranded_counts[base_to_idx_stranded[alt_base.lower()]]) / float(sum(stranded_counts[4:]))
    dp4 = [ stranded_counts[base_to_idx_stranded[ref_base]],
            stranded_counts[base_to_idx_stranded[ref_base.lower()]],
            stranded_counts[base_to_idx_stranded[alt_base]],
            stranded_counts[base_to_idx_stranded[alt_base.lower()]]]
    
    return (dp,af,faf,raf,dp4)

def parseIndelPileup(fields,ref_base,alt_base):
    counts = [0,0,0,0,0,0,0,0,0] # indel ref match, indel fwd ref match, indel rev ref match, indel alt match, indel fwd alt match, indel rev alt match, other, other fwd, other rev
    ref_base2 = fields[2]

    carrot_flag = False
    ins_flag = False
    ins_str = ""
    ins_len = 0
    del_flag = False
    del_str = ""
    del_len = 0
    first_iter = True
    forward_flag = False
    last_seq = ""
    last_seq_code = 'b'
    for base in fields[4]:
        if ins_flag:
            if base.isdigit():
                ins_str += base
            else:
                ins_len = int(ins_str)
                ins_flag = False
        if del_flag:
            if base.isdigit():
                del_str += base
            else:
                del_len = int(del_str)
                del_flag = False
        if ins_len > 0:
            last_seq += base
            last_seq_code = 'i'
            ins_len -= 1
            continue
        if del_len > 0:
            last_seq += base
            last_seq_code = 'd'
            del_len -= 1
            continue
        if carrot_flag:
            carrot_flag = False
            continue
        if base == '.' or base == ','\
            or base == 'A' or base == 'a'\
            or base == 'C' or base == 'c'\
            or base == 'G' or base == 'g'\
            or base == 'T' or base == 't'\
            or base == 'N' or base == 'n'\
            or base == '>' or base == '<'\
            or base == '*' or base == '#':
            if first_iter:
                first_iter = False
            else:
                if last_seq_code == 'i':
                    if last_seq.upper() == alt_base.upper():
                        counts[3] += 1
                        if forward_flag:
                            counts[4] += 1
                        else:
                            counts[5] += 1
                    else:
                        counts[6] += 1
                        if forward_flag:
                            counts[7] += 1
                        else:
                            counts[8] += 1
                elif last_seq_code == 'd':
                    if last_seq.upper() == ref_base.upper():
                        counts[3] += 1
                        if forward_flag:
                            counts[4] += 1
                        else:
                            counts[5] += 1
                    else:
                        counts[6] += 1
                        if forward_flag:
                            counts[7] += 1
                        else:
                            counts[8] += 1
                elif last_seq_code == 'b':
                    if last_seq.upper() == ref_base.upper():
                        counts[0] += 1
                        if forward_flag:
                            counts[1] += 1
                        else:
                            counts[2] += 1
                    elif last_seq.upper() == alt_base.upper():
                        counts[3] += 1
                        if forward_flag:
                            counts[4] += 1
                        else:
                            counts[5] += 1
                    else:
                        counts[6] += 1
                        if forward_flag:
                            counts[7] += 1
                        else:
                            counts[8] += 1
            if base == '.':
                last_seq = ref_base2
                forward_flag = True
                last_seq_code = 'b'
            elif base == ',':
                last_seq = ref_base2
                forward_flag = False
                last_seq_code = 'b'
            elif base == '>' or base == '<' or base == '*' or base == '#':
                continue
            else:
                forward_flag = base.isupper()
                last_seq = base.upper()
                last_seq_code = 'b'
        elif base == '+':
            ins_flag = True
            ins_str = ""
        elif base == '-':
            del_flag = True
            del_str = ""
        elif base == '^':
            carrot_flag = True
        elif base == '$':
            continue
        if first_iter:
            first_iter = False

    if last_seq_code == 'i':
        if last_seq.upper() == alt_base.upper():
            counts[3] += 1
            if forward_flag:
                counts[4] += 1
            else:
                counts[5] += 1
        else:
            counts[6] += 1
            if forward_flag:
                counts[7] += 1
            else:
                counts[8] += 1
    elif last_seq_code == 'd':
        if last_seq.upper() == ref_base.upper():
            counts[3] += 1
            if forward_flag:
                counts[4] += 1
            else:
                counts[5] += 1
        else:
            counts[6] += 1
            if forward_flag:
                counts[7] += 1
            else:
                counts[8] += 1
    elif last_seq_code == 'b':
        if last_seq.upper() == ref_base.upper():
            counts[0] += 1
            if forward_flag:
                counts[1] += 1
            else:
                counts[2] += 1
        elif last_seq.upper() == alt_base.upper():
            counts[3] += 1
            if forward_flag:
                counts[4] += 1
            else:
                counts[5] += 1
        else:
            counts[6] += 1
            if forward_flag:
                counts[7] += 1
            else:
                counts[8] += 1
    dp = int(fields[3])
    af = float(counts[3]) / float(sum([counts[0],counts[3],counts[6]]))
    if sum([counts[1],counts[4],counts[7]]) == 0:
        faf = float("nan")
    else:
        faf = float(counts[4]) / float(sum([counts[1],counts[4],counts[7]]))
    if sum([counts[2],counts[5],counts[8]]) == 0:
        raf = float("nan")
    else:
        raf = float(counts[5]) / float(sum([counts[2],counts[5],counts[8]]))
    dp4 = [counts[1],counts[2],counts[4],counts[5]]
    return (dp,af,faf,raf,dp4)

def annotateVCF(in_vcf_filepath,in_mpileup_filepath,out_vcf_filepath):
    in_vcf = open(in_vcf_filepath, 'r')
    in_mpileup = open(in_mpileup_filepath,'r')
    out_vcf = open(out_vcf_filepath,'w')

    #First pass parsing of input vcf, output headerlines + new headerlines, add VCF sites we care about to to_examine (limits memory usage for sites that don't need annotation)
    to_examine = {}
    for line in in_vcf:
        if line[0:2] == "##":
            out_vcf.write(line)
        elif line[0] == "#":
            out_vcf.write("##annotateVCFVersion=0.1\n")
            out_vcf.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">\n")
            out_vcf.write("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n")
            out_vcf.write("##INFO=<ID=FAF,Number=1,Type=Float,Description=\"Forward Allele Frequency\">\n")
            out_vcf.write("##INFO=<ID=RAF,Number=1,Type=Float,Description=\"Reverse Allele Frequency\">\n")
            out_vcf.write("##INFO=<ID=SB,Number=1,Type=Integer,Description=\"Phred-scaled strand bias at this position\">\n")
            out_vcf.write("##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">\n")
            out_vcf.write(line)
        else:
            fields = line.strip().split()
            if fields[0] in to_examine:
                to_examine[fields[0]][int(fields[1])] = (fields[3],fields[4])
            else:
                to_examine[fields[0]] = {int(fields[1]) : (fields[3],fields[4])}
    in_vcf.close()
    data = {}

    #Populate data dictionary, which relates chromosome and position to the following:
    # depth of coverage
    # allele frequency
    # forward strand allele frequency
    # reverse strand allele frequency
    # dp4 - depth of coverage of ref allele fwd strand, DOC of ref allele rev strand, DOC of alt allele fwd strand, DOC of alt allele rev strand
    for line in in_mpileup:
        fields = line.strip().split()
        if fields[0] not in to_examine:
            continue
        if int(fields[1]) not in to_examine[fields[0]]: #
            continue
        (ref_base,alt_base) = to_examine[fields[0]][int(fields[1])]
        if len(ref_base.split(',')) > 1: # Can't handle multiple ref alleles
            continue
        if len(alt_base.split(',')) > 1: # Can't handle multiple alt alleles
            continue
        if len(ref_base) > 1 or len(alt_base) > 1:
            if len(ref_base) > 1 and len(alt_base) > 1: # Can't handle complex indels
                continue
            data[(fields[0],int(fields[1]))] = parseIndelPileup(fields,ref_base,alt_base)
        if len(ref_base) == 1 and len(alt_base) == 1:
            data[(fields[0],int(fields[1]))] = parseSimpleSNPpileup(fields,ref_base,alt_base)
        
    in_mpileup.close()
    #Reopen vcf, this time, skip header, annotate all the sites for which there is an entry in data dictionary
    #(Sites without entries have either multiple ref or alt bases, or have complex indels. Not supported (for now), and not reported as a result)
    in_vcf = open(in_vcf_filepath,'r')
    for line in in_vcf:
        if line[0] == '#':
            continue
        fields = line.strip().split('\t')
        if (fields[0],int(fields[1])) not in data:
            continue
        (dp,af,faf,raf,dp4) = data[(fields[0],int(fields[1]))]
        dp2x2 = [[dp4[0],dp4[1]],[dp4[2],dp4[3]]]
        _, p_val = fisher_exact(dp2x2)
        sb = pval_to_phredqual(p_val)
        if fields[7] == "":
            info = []
        else:
            info = fields[7].split(';')
        info.append("DP=%d" %(dp))
        info.append("AF=%.6f" %(af))
        if isnan(faf):
            info.append("FAF=NaN")
        else:
            info.append("FAF=%.6f" %(faf))
        if isnan(raf):
            info.append("RAF=NaN")
        else:
            info.append("RAF=%.6f" %(raf))
        info.append("SB=%d" %(sb))
        info.append("DP4=%s"%(','.join([str(x) for x in dp4])))
        new_info = ';'.join(info)
        fields[7] = new_info
        out_vcf.write("%s\n" %("\t".join(fields)))
    in_vcf.close()
    out_vcf.close()

if __name__ == "__main__":
    annotateVCF(sys.argv[1],sys.argv[2],sys.argv[3]) #