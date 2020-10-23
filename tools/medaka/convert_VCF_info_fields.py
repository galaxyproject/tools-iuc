#!/usr/bin/env python3

# Takes in VCF file annotated with medaka tools annotate and converts
#
# Usage statement:
# python convert_VCF_info_fields.py in_vcf.vcf out_vcf.vcf

# 10/21/2020 - Nathan P. Roach, natproach@gmail.com

import sys
from math import log10

from scipy.stats import fisher_exact


def pval_to_phredqual(pval):
    return int(round(-10. * log10(pval)))


def parseInfoField(info):
    info_fields = info.split(';')
    info_dict = {}
    info_list = []  # Want to keep fields in order they appear in original VCF
    for info_field in info_fields:
        code, val = info_field.split('=')
        info_dict[code] = val
        info_list.append(code)
    return info_dict, info_list


def annotateVCF(in_vcf_filepath, out_vcf_filepath):
    in_vcf = open(in_vcf_filepath, 'r')
    out_vcf = open(out_vcf_filepath, 'w')
    to_skip = set()
    for i, line in enumerate(in_vcf):
        if i == 1:
            out_vcf.write("##convert_VCF_info_fields=0.1\n")
        if line[0:2] == "##":
            if "Number=." in line:
                assert line[0:11] == "##INFO=<ID="
                line = line[11:]
                id_ = line.split(',')[0]
                to_skip.add(id_)
                continue
            out_vcf.write(line)
        elif line[0] == "#":
            # out_vcf.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">\n") # Already annotated by medaka tools annotate
            out_vcf.write("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n")
            out_vcf.write("##INFO=<ID=FAF,Number=1,Type=Float,Description=\"Forward Allele Frequency\">\n")
            out_vcf.write("##INFO=<ID=RAF,Number=1,Type=Float,Description=\"Reverse Allele Frequency\">\n")
            out_vcf.write("##INFO=<ID=SB,Number=1,Type=Integer,Description=\"Phred-scaled strand bias at this position\">\n")
            out_vcf.write("##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">\n")
            out_vcf.write("##INFO=<ID=AS,Number=4,Type=Integer,Description=\"Total alignment score to ref and alt allele of spanning reads by strand (ref fwd, ref rev, alt fwd, alt rev) aligned with parasail match 5, mismatch -4, open 5, extend 3\">\n")
            out_vcf.write(line)
        else:
            fields = line.strip().split('\t')
            info_dict, info_list = parseInfoField(fields[7])
            sr_list = [int(x) for x in info_dict["SR"].split(',')]
            sc_list = [int(x) for x in info_dict["SC"].split(',')]
            if len(sr_list) == len(sc_list):
                print(len(sr_list))
                variant_count = int((len(sr_list) / 2) - 1)
                dp = int(info_dict["DP"])
                dpf = int(info_dict["DPS"].split(',')[0])
                dpr = int(info_dict["DPS"].split(',')[1])
                for i in range(variant_count):
                    ref_fwd = 0
                    ref_rev = 1
                    var_fwd = 2 * (i + 1)
                    var_rev = 2 * (i + 1) + 1

                    dp4 = [sr_list[ref_fwd], sr_list[ref_rev], sr_list[var_fwd], sr_list[var_rev]]
                    dp2x2 = [[dp4[0], dp4[1]], [dp4[2], dp4[3]]]
                    _, p_val = fisher_exact(dp2x2)
                    sb = pval_to_phredqual(p_val)

                    as_ = [sc_list[ref_fwd], sc_list[ref_rev], sc_list[var_fwd], sc_list[var_rev]]

                    info = []
                    for code in info_list:
                        val = info_dict[code]
                        info.append("%s=%s" % (code, val))

                    if dp == 0:
                        info.append("AF=NaN")
                    else:
                        af = float(dp4[2] + dp4[3]) / float(dp)
                        info.append("AF=%.6f" % (af))
                    if dpf == 0:
                        info.append("FAF=NaN")
                    else:
                        faf = float(dp4[2]) / float(dpf)
                        info.append("FAF=%.6f" % (faf))
                    if dpr == 0:
                        info.append("RAF=NaN")
                    else:
                        raf = float(dp4[3]) / float(dpr)
                        info.append("RAF=%.6f" % (raf))
                    info.append("SB=%d" % (sb))
                    info.append("DP4=%s" % (','.join([str(x) for x in dp4])))
                    info.append("AS=%s" % (','.join([str(x) for x in as_])))
                    new_info = ';'.join(info)
                    fields[7] = new_info
                    out_vcf.write("%s\n" % ("\t".join(fields)))
            else:
                print("WARNING - SR and SC are different lengths, skipping variant")
                continue
    in_vcf.close()
    out_vcf.close()


if __name__ == "__main__":
    annotateVCF(sys.argv[1], sys.argv[2])
