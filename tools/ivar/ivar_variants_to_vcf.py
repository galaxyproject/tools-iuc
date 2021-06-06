#!/usr/bin/env python
import os
import sys
import re
import errno
import argparse


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
    parser.add_argument(
        "-af",
        "--allele_freq_thresh",
        type=float,
        dest="ALLELE_FREQ_THRESH",
        default=0,
        help="Only output variants where allele frequency greater than this number (default: 0).",
    )

    return parser.parse_args(args)


def make_dir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def ivar_variants_to_vcf(FileIn, FileOut, passOnly=False, minAF=0):
    filename = os.path.splitext(FileIn)[0]
    header = (
        "##fileformat=VCFv4.2\n"
        "##source=iVar\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'
        '##FILTER=<ID=PASS,Description="Result of p-value <= 0.05">\n'
        '##FILTER=<ID=FAIL,Description="Result of p-value > 0.05">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=REF_DP,Number=1,Type=Integer,Description="Depth of reference base">\n'
        '##FORMAT=<ID=REF_RV,Number=1,Type=Integer,Description="Depth of reference base on reverse reads">\n'
        '##FORMAT=<ID=REF_QUAL,Number=1,Type=Integer,Description="Mean quality of reference base">\n'
        '##FORMAT=<ID=ALT_DP,Number=1,Type=Integer,Description="Depth of alternate base">\n'
        '##FORMAT=<ID=ALT_RV,Number=1,Type=Integer,Description="Deapth of alternate base on reverse reads">\n'
        '##FORMAT=<ID=ALT_QUAL,Number=1,Type=String,Description="Mean quality of alternate base">\n'
        '##FORMAT=<ID=ALT_FREQ,Number=1,Type=String,Description="Frequency of alternate base">\n'
    )
    header += (
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + filename + "\n"
    )

    varList = []
    varCountDict = {"SNP": 0, "INS": 0, "DEL": 0}
    OutDir = os.path.dirname(FileOut)
    make_dir(OutDir)
    fout = open(FileOut, "w")
    fout.write(header)
    with open(FileIn) as f:
        for line in f:
            if not re.match("REGION", line):
                line = re.split("\t", line)
                CHROM = line[0]
                POS = line[1]
                ID = "."
                REF = line[2]
                ALT = line[3]
                var_type = "SNP"
                if ALT[0] == "+":
                    ALT = REF + ALT[1:]
                    var_type = "INS"
                elif ALT[0] == "-":
                    REF += ALT[1:]
                    ALT = line[2]
                    var_type = "DEL"
                QUAL = "."
                pass_test = line[13]
                if pass_test == "TRUE":
                    FILTER = "PASS"
                else:
                    FILTER = "FAIL"
                INFO = "DP=" + line[11]
                FORMAT = "GT:REF_DP:REF_RV:REF_QUAL:ALT_DP:ALT_RV:ALT_QUAL:ALT_FREQ"
                SAMPLE = (
                    "1:"
                    + line[4]
                    + ":"
                    + line[5]
                    + ":"
                    + line[6]
                    + ":"
                    + line[7]
                    + ":"
                    + line[8]
                    + ":"
                    + line[9]
                    + ":"
                    + line[10]
                )
                oline = (
                    CHROM
                    + "\t"
                    + POS
                    + "\t"
                    + ID
                    + "\t"
                    + REF
                    + "\t"
                    + ALT
                    + "\t"
                    + QUAL
                    + "\t"
                    + FILTER
                    + "\t"
                    + INFO
                    + "\t"
                    + FORMAT
                    + "\t"
                    + SAMPLE
                    + "\n"
                )
                writeLine = True
                if passOnly and FILTER != "PASS":
                    writeLine = False
                if float(line[10]) < minAF:
                    writeLine = False
                if (CHROM, POS, REF, ALT) in varList:
                    writeLine = False
                else:
                    varList.append((CHROM, POS, REF, ALT))
                if writeLine:
                    varCountDict[var_type] += 1
                    fout.write(oline)
    fout.close()

    # Print variant counts to pass to MultiQC
    varCountList = [(k, str(v)) for k, v in sorted(varCountDict.items())]
    print("\t".join(["sample"] + [x[0] for x in varCountList]))
    print("\t".join([filename] + [x[1] for x in varCountList]))


def main(args=None):
    args = parse_args(args)
    ivar_variants_to_vcf(
        args.FILE_IN, args.FILE_OUT, args.PASS_ONLY, args.ALLELE_FREQ_THRESH
    )


if __name__ == "__main__":
    sys.exit(main())
