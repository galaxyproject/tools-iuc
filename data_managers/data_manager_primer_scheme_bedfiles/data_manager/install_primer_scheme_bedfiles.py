#!/usr/bin/env python

from __future__ import division, print_function

import argparse
try:
    from io import StringIO
except ImportError:
    from StringIO import StringIO
import json
import os
import os.path
import re
import sys

import requests

DATA_TABLE_NAME = "primer_scheme_bedfiles"


def convert_and_write_bed(input_file, bed_output_filename, scheme_name, force_string=True):
    with open(bed_output_filename, "w") as bed_output_file:
        for line in input_file:
            fields = line.strip().split("\t")
            if "Midnight" in scheme_name:
                # Midnight primers are distributed in a tabular file, not a BED file
                if line.startswith("Primer Name"):
                    continue
                if len(fields) != 8:
                    exit("Unexpected format in Midnight primer file: {}".format(line.rstrip()))
                (primer_name, _, pool, _, _, _, start, end) = fields
                strand = '+' if primer_name.endswith('LEFT') else '-'
                if strand == '-':
                    start, end = end, start
                fields = ["MN908947.3", start, end, primer_name, pool, strand]
            else:
                if len(fields) < 5:
                    # too short to encode the "ARTIC style BED" format
                    exit("invalid format in BED file: {}".format(line.rstrip()))
            # 'BED' format used by ARTIC pipeline uses
            # chrom  start  end  primer_name  pool_name
            # see this: https://github.com/artic-network/fieldbioinformatics/blob/master/artic/vcftagprimersites.py#L76
            # for ARTIC minion and
            # this: https://github.com/andersen-lab/ivar/blob/master/src/primer_bed.cpp#L125
            # for ivar trim (ivar trim treats the file as BED following the standard but also allows the ARTIC format)
            try:
                float(fields[4])
            except ValueError:
                # this is a string, we can leave it as is
                pass
            else:
                # ensure that it is forced to be a string
                fields[4] = '_{0}'.format(fields[4])
            print('\t'.join(fields), file=bed_output_file)


def fetch_primers(output_directory, primers):
    primer_sets = {
        "SARS-CoV-2-ARTICv1": "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V1/nCoV-2019.bed",
        "SARS-CoV-2-ARTICv2": "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V2/nCoV-2019.bed",
        "SARS-CoV-2-ARTICv3": "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V3/nCoV-2019.bed",
        "SARS-CoV-2-ARTICv4": "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V4/SARS-CoV-2.scheme.bed",
        "VarSkip-V1a": "https://raw.githubusercontent.com/nebiolabs/VarSkip/main/schemes/NEB_VarSkip/V1a/NEB_VarSkip.scheme.bed",
        "Midnight-v1": "https://zenodo.org/record/3897530/files/SARS-CoV-2_primer_sets_RBK004_nanopore_sequencing.tab?download=1"
    }

    data = []
    for name, url in primer_sets.items():
        if name not in primers:
            continue
        response = requests.get(url)
        if response.status_code != 200:
            print(
                "Error: download of",
                url,
                "failed with code",
                response.status_code,
                file=sys.stderr,
            )
            exit(response.status_code)
        bed_output_filename = os.path.join(output_directory, name + ".bed")
        convert_and_write_bed(StringIO(response.text), bed_output_filename, name)
        if 'ARTIC' in name:
            # split the vX from the rest of the name in ARTIC primer set description
            description = name[:-2] + " " + name[-2:] + " primer set"
        else:
            description = name + " primer set"
        data.append(dict(value=name, path=bed_output_filename, description=description))
    return data


def install_primer_file(
    output_directory, input_filename, scheme_name, primer_description
):
    name = re.sub(r"[^\w-]", "", str(scheme_name).replace(" ", "_"))
    output_filename = os.path.join(output_directory, name + ".bed")
    with open(input_filename) as input_file:
        convert_and_write_bed(input_file, output_filename, scheme_name)
    data = [dict(value=name, description=primer_description, path=output_filename)]
    return data


class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(","))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch ARTIC, VarSkip and Midnight SARS-CoV-2 primer files for Galaxy/IRIDA use"
    )
    parser.add_argument(
        "--output_directory", default="tmp", help="Directory to write output to"
    )
    primer_file = parser.add_argument_group()
    primer_file.add_argument(
        "--primer_file", help="BED format file containing primer scheme"
    )
    primer_file.add_argument(
        "--primer_name",
        help="Name of primer scheme (one word). Required if --primer_file is used",
    )
    primer_file.add_argument(
        "--primer_description",
        help="Description of primer scheme. Required if --primer_file is used",
    )
    artic = parser.add_argument_group()
    artic.add_argument(
        "--artic_primers",
        action=SplitArgs,
        help="Comma separated list of primers to fetch",
    )
    parser.add_argument(
        "galaxy_datamanager_filename",
        help="Galaxy JSON format file describing data manager inputs",
    )
    args = parser.parse_args()

    if args.artic_primers is None and args.primer_file is None:
        print(
            "One of --artic_primers or --primer_file + --primer_name + --primer_description is required.",
            file=sys.stderr,
        )
        exit(1)
    elif args.primer_file is not None and (
        args.primer_name is None or args.primer_description is None
    ):
        print(
            "If --primer_file is used --primer_name and --primer_description is also required",
            file=sys.stderr,
        )
        exit(1)
    elif args.primer_file is not None and args.artic_primers is not None:
        print(
            "Only one of --artic_primers or --primer_file + --primer_name + --primer_description can be chosen"
        )
        exit(1)

    with open(args.galaxy_datamanager_filename) as fh:
        config = json.load(fh)
    output_directory = config.get("output_data", [{}])[0].get("extra_files_path", None)
    if output_directory is None:
        output_directory = args.output_directory

    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    data_manager_dict = {}
    data_manager_dict["data_tables"] = config.get("data_tables", {})
    data_manager_dict["data_tables"][DATA_TABLE_NAME] = []

    if args.artic_primers:
        data = fetch_primers(output_directory, args.artic_primers)
    else:
        data = install_primer_file(
            output_directory,
            args.primer_file,
            args.primer_name,
            args.primer_description,
        )

    data_manager_dict["data_tables"][DATA_TABLE_NAME].extend(data)
    with open(args.galaxy_datamanager_filename, "w") as fh:
        json.dump(data_manager_dict, fh, sort_keys=True)
