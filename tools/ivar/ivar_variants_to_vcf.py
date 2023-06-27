#!/usr/bin/env python
import argparse
import errno
import os
import re
import sys


def parse_args(args=None):
    Description = "Convert iVar variants tsv file to vcf format."
    Epilog = """Example usage: python ivar_variants_to_vcf.py <FILE_IN> <FILE_OUT>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input tsv file.")
    parser.add_argument("FILE_OUT", help="Full path to output vcf file.")
    parser.add_argument(
        "-po",
        "--pass_only",
        dest="PASS_ONLY",
        help="Only output variants that PASS all filters.",
        action="store_true",
    )

    return parser.parse_args(args)


def make_dir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def info_line(info_keys, kv):
    info_strings = []
    for key in info_keys:
        if key not in kv:
            raise KeyError(
                'Expected key {} missing from INFO field key value pairs'
                .format(key)
            )
        if kv[key] is False:
            # a FLAG element, which should not be set
            continue
        if kv[key] is True:
            # a FLAG element => write the key only
            info_strings.append(key)
        else:
            info_strings.append('{}={}'.format(key, kv[key]))
    return ';'.join(info_strings)


def ivar_variants_to_vcf(FileIn, FileOut, passOnly=False):
    filename = os.path.splitext(FileIn)[0]
    header = (
        "##fileformat=VCFv4.2\n"
        "##source=iVar\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'
        '##INFO=<ID=REF_DP,Number=1,Type=Integer,Description="Depth of reference base">\n'
        '##INFO=<ID=REF_RV,Number=1,Type=Integer,Description="Depth of reference base on reverse reads">\n'
        '##INFO=<ID=REF_QUAL,Number=1,Type=Integer,Description="Mean quality of reference base">\n'
        '##INFO=<ID=ALT_DP,Number=1,Type=Integer,Description="Depth of alternate base">\n'
        '##INFO=<ID=ALT_RV,Number=1,Type=Integer,Description="Deapth of alternate base on reverse reads">\n'
        '##INFO=<ID=ALT_QUAL,Number=1,Type=Integer,Description="Mean quality of alternate base">\n'
        '##INFO=<ID=AF,Number=1,Type=Float,Description="Frequency of alternate base">\n'
        '##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">\n'
        '##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">'
        '##FILTER=<ID=PASS,Description="Result of p-value <= 0.05">\n'
        '##FILTER=<ID=FAIL,Description="Result of p-value > 0.05">\n'
    )
    info_keys = [
        re.match(r'##INFO=<ID=([^,]+),', line).group(1)
        for line in header.splitlines()
        if line.startswith('##INFO=')
    ]
    header += (
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )

    vars_seen = set()
    varCountDict = {"SNP": 0, "INS": 0, "DEL": 0}
    OutDir = os.path.dirname(FileOut)
    make_dir(OutDir)
    with open(FileIn) as f, open(FileOut, "w") as fout:
        fout.write(header)
        for line in f:
            if line.startswith("REGION"):
                continue

            # fields: 
            # 0 REGION
            # 1 POS
            # 2 REF
            # 3 ALT
            # 4 REF_DP
            # 5 REF_RV
            # 6 REF_QUAL
            # 7 ALT_DP
            # 8 ALT_RV
            # 9 ALT_QUAL
            # 10 ALT_FREQ
            # 11 TOTAL_DP
            # 12 PVAL
            # 13 PASS
            # 14 GFF_FEATURE
            # 15 REF_CODON
            # 16 REF_AA
            # 17 ALT_CODON
            # 18 ALT_AA
            # 19 POS_AA
            line = line.split("\t")
            CHROM = line[0]
            POS = line[1]
            ID = "."
            REF = line[2]
            ALT = line[3]
            if ALT[0] == "+":
                ALT = REF + ALT[1:]
                var_type = "INS"
            elif ALT[0] == "-":
                REF += ALT[1:]
                ALT = line[2]
                var_type = "DEL"
            else:
                var_type = "SNP"
            QUAL = "."
            pass_test = line[13]
            if pass_test == "TRUE":
                FILTER = "PASS"
            else:
                FILTER = "FAIL"

            if (passOnly and FILTER != "PASS"):
                continue
            var = (CHROM, POS, REF, ALT)
            if var in vars_seen:
                continue

            ref_dp = int(line[4])
            ref_dp_rev = int(line[5])
            ref_dp_fwd = ref_dp - ref_dp_rev
            
            alt_dp = int(line[7])
            alt_dp_rev = int(line[8])
            alt_dp_fwd = alt_dp - alt_dp_rev

            dp4 = f'{ref_dp_fwd},{ref_dp_rev},{alt_dp_fwd},{alt_dp_rev}'
            info_elements = {
                'DP': line[11],
                'REF_DP': ref_dp,
                'REF_RV': ref_dp_rev,
                'REF_QUAL': line[6],
                'ALT_DP': alt_dp,
                'ALT_RV': alt_dp_rev,
                'ALT_QUAL': line[9],
                'AF': line[10],
                'DP4': dp4
            }
            if var_type in ['INS', 'DEL']:
                # add INDEL FLAG
                info_elements['INDEL'] = True
            else:
                info_elements['INDEL'] = False
            INFO = info_line(info_keys, info_elements)

            vars_seen.add(var)
            varCountDict[var_type] += 1
            fout.write(
                '\t'.join(
                    [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO]
                ) + '\n'
            )

    # Print variant counts to pass to MultiQC
    varCountList = [(k, str(v)) for k, v in sorted(varCountDict.items())]
    print("\t".join(["sample"] + [x[0] for x in varCountList]))
    print("\t".join([filename] + [x[1] for x in varCountList]))


def main(args=None):
    args = parse_args(args)
    ivar_variants_to_vcf(
        args.FILE_IN, args.FILE_OUT, args.PASS_ONLY
    )


if __name__ == "__main__":
    sys.exit(main())
