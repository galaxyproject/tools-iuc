#!/usr/bin/env python3

# Takes in VCF file annotated with medaka tools annotate and converts
#
# Usage statement:
# python convert_VCF_info_fields.py in_vcf.vcf out_vcf.vcf

# 10/21/2020 - Nathan P. Roach, natproach@gmail.com

import sys
from collections import OrderedDict
from math import log10

from scipy.stats import fisher_exact


def pval_to_phredqual(pval):
    try:
        ret = round(-10 * log10(pval))
    except ValueError:
        ret = 2147483647  # transform pval of 0.0 to max signed 32 bit int
    return ret


def parseInfoField(info):
    info_fields = info.split(';')
    info_dict = OrderedDict()
    for info_field in info_fields:
        code, val = info_field.split('=')
        info_dict[code] = val
    return info_dict


def annotateVCF(in_vcf_filepath, out_vcf_filepath):
    in_vcf = open(in_vcf_filepath, 'r')
    out_vcf = open(out_vcf_filepath, 'w')
    to_skip = set(['SC', 'SR'])
    for i, line in enumerate(in_vcf):
        if i == 1:
            out_vcf.write("##convert_VCF_info_fields=0.2\n")
        if line[0:2] == "##":
            if line[0:11] == "##INFO=<ID=":
                id_ = line[11:].split(',')[0]
                if id_ in to_skip:
                    continue
            out_vcf.write(line)
        elif line[0] == "#":
            out_vcf.write('##INFO=<ID=DPSPS,Number=2,Type=Integer,Description="Spanning Reads Allele Frequency By Strand">\n')
            out_vcf.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Spanning Reads Allele Frequency">\n')
            out_vcf.write('##INFO=<ID=FAF,Number=1,Type=Float,Description="Forward Spanning Reads Allele Frequency">\n')
            out_vcf.write('##INFO=<ID=RAF,Number=1,Type=Float,Description="Reverse Spanning Reads Allele Frequency">\n')
            out_vcf.write('##INFO=<ID=SB,Number=1,Type=Integer,Description="Phred-scaled strand bias of spanning reads at this position">\n')
            out_vcf.write('##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases in spanning reads">\n')
            out_vcf.write('##INFO=<ID=AS,Number=4,Type=Integer,Description="Total alignment score to ref and alt allele of spanning reads by strand (ref fwd, ref rev, alt fwd, alt rev) aligned with parasail match 5, mismatch -4, open 5, extend 3">\n')
            out_vcf.write(line)
        else:
            fields = line.split('\t')
            info_dict = parseInfoField(fields[7])
            sr_list = [int(x) for x in info_dict["SR"].split(',')]
            sc_list = [int(x) for x in info_dict["SC"].split(',')]
            if len(sr_list) == len(sc_list):
                variant_list = fields[4].split(',')
                dpsp = int(info_dict["DPSP"])
                ref_fwd, ref_rev = 0, 1
                dpspf, dpspr = (int(x) for x in info_dict["AR"].split(','))
                for i in range(0, len(sr_list), 2):
                    dpspf += sr_list[i]
                    dpspr += sr_list[i + 1]
                for j, i in enumerate(range(2, len(sr_list), 2)):
                    dp4 = (sr_list[ref_fwd], sr_list[ref_rev], sr_list[i], sr_list[i + 1])
                    dp2x2 = [[dp4[0], dp4[1]], [dp4[2], dp4[3]]]
                    _, p_val = fisher_exact(dp2x2)
                    sb = pval_to_phredqual(p_val)

                    as_ = (sc_list[ref_fwd], sc_list[ref_rev], sc_list[i], sc_list[i + 1])

                    info = []
                    for code in info_dict:
                        if code in to_skip:
                            continue
                        val = info_dict[code]
                        info.append("%s=%s" % (code, val))

                    info.append("DPSPS=%d,%d" % (dpspf, dpspr))

                    if dpsp == 0:
                        info.append("AF=NaN")
                    else:
                        af = (dp4[2] + dp4[3]) / dpsp
                        info.append("AF=%.6f" % (af))
                    if dpspf == 0:
                        info.append("FAF=NaN")
                    else:
                        faf = dp4[2] / dpspf
                        info.append("FAF=%.6f" % (faf))
                    if dpspr == 0:
                        info.append("RAF=NaN")
                    else:
                        raf = dp4[3] / dpspr
                        info.append("RAF=%.6f" % (raf))
                    info.append("SB=%d" % (sb))
                    info.append("DP4=%d,%d,%d,%d" % (dp4))
                    info.append("AS=%d,%d,%d,%d" % (as_))
                    new_info = ';'.join(info)
                    fields[4] = variant_list[j]
                    fields[7] = new_info
                    out_vcf.write("%s" % ("\t".join(fields)))
            else:
                print("WARNING - SR and SC are different lengths, skipping variant")
                print(line.strip())  # Print the line for debugging purposes
    in_vcf.close()
    out_vcf.close()


if __name__ == "__main__":
    annotateVCF(sys.argv[1], sys.argv[2])
