import argparse
from collections import OrderedDict

from pyfaidx import Fasta

parser = argparse.ArgumentParser(description='Generate report tables')
parser.add_argument("--reference",
                    required=True,
                    help="Filepath to reference FASTA file")
parser.add_argument("--in-vcf",
                    required=True,
                    help="Filepath to vcf file to be analyzed")
parser.add_argument("--out-vcf",
                    required=True,
                    help="Filepath to vcf file to be output")


def parseInfoField(info):
    info_fields = info.split(';')
    info_dict = OrderedDict()
    for info_field in info_fields:
        info_field_fields = info_field.split('=')
        if len(info_field_fields) == 2:
            code, val = info_field_fields
        else:
            code = info_field_fields[0]
            val = ""
        info_dict[code] = val
    return info_dict


args = parser.parse_args()
ref_path = args.reference
reference = Fasta(ref_path, sequence_always_upper=True, read_ahead=1000)
in_vcf_path = args.in_vcf
in_vcf = open(in_vcf_path)
out_vcf_path = args.out_vcf
out_vcf = open(out_vcf_path, 'w')

first_info = True
for line in in_vcf:
    if line[0] == '#':
        if len(line) >= 2:
            if line[0:2] == "##":
                if len(line) >= 6:
                    if first_info and line[0:6] == "##INFO":
                        out_vcf.write("##INFO=<ID=HRUN,Number=1,Type=Integer,Description=\"Homopolymer length to the right of report indel position\">\n")
                        first_info = False
            elif first_info:  # There were no INFO lines; add one before fields header
                out_vcf.write("##INFO=<ID=HRUN,Number=1,Type=Integer,Description=\"Homopolymer length to the right of report indel position\">\n")
        out_vcf.write(line)
    else:
        fields = line.split('\t')
        chrom = fields[0]
        pos = int(fields[1]) - 1
        ref = fields[3]
        alt = fields[4]
        if len(ref) != len(alt):  # Indel
            info_dict = parseInfoField(fields[7])

            window = 50
            hrun = 1
            start = pos + 2
            end = start + window
            base = reference[chrom][pos + 1]
            seq_len = len(reference[chrom])
            for i in range(start, len(reference)):
                base2 = reference[chrom][i]
                if base == base2:
                    hrun += 1
                else:
                    break
            # Extend to left in case not left aligned
            for i in range(pos, -1, -1):
                if reference[chrom][i] == base:
                    hrun += 1
                else:
                    break
            info = ["HRUN=%d" % (hrun)]
            for code in info_dict:
                val = info_dict[code]
                if val == "":
                    info.append("%s" % (code))
                else:
                    info.append("%s=%s" % (code, val))
            new_info = ';'.join(info)
            fields[7] = new_info
            out_vcf.write("%s" % ("\t".join(fields)))
        else:
            out_vcf.write(line)
in_vcf.close()
out_vcf.close()
