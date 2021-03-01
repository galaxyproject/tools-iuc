#!/usr/bin/env python
import argparse
from collections import OrderedDict

import vcf
from pyfaidx import Fasta
from vcf.parser import _Info as VcfInfo


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
in_vcf = vcf.Reader(open(in_vcf_path))
in_vcf.infos['HRUN'] = VcfInfo(
    'HRUN', 1, 'Integer',
    'Homopolymer length to the right of report indel position',
    "get_hrun", "1.0")
out_vcf_path = args.out_vcf
out_vcf = vcf.Writer(open(out_vcf_path, 'w'), in_vcf)
for record in in_vcf:
    chrom = record.CHROM
    pos = record.POS - 1
    ref = record.REF
    calc_hrun = False
    for alt in record.ALT:
        if len(ref) != len(alt):
            calc_hrun = True
    if calc_hrun:
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
        record.add_info('HRUN', [hrun])
    out_vcf.write_record(record)
# in_vcf.close()
out_vcf.close()
