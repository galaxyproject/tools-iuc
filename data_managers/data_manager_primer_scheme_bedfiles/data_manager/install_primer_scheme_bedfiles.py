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


def write_artic_style_bed(input_file, bed_output_filename):
    with open(bed_output_filename, "w") as bed_output_file:
        for line in input_file:
            fields = line.split("\t")
            if len(fields) < 6:
                # too short to encode the strand format
                exit("invalid format in BED file: {}".format(line.rstrip()))
            try:
                # try and parse field 5 as a number
                score = float(fields[4])
            except ValueError:
                # Alright, this is an ARTIC-style bed,
                # which is actually against the specs, but required by the
                # ARTIC pipeline.
                pass
            else:
                # This is a regular bed with numbers in the score column.
                # We need to "fix" it for the ARTIC pipeline.
                fields[4] = '_{0}'.format(score)
            bed_output_file.write("\t".join(fields))


def fetch_artic_primers(output_directory, primers):
    primer_sets = {
        "SARS-CoV-2-ARTICv1": "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V1/nCoV-2019.bed",
        "SARS-CoV-2-ARTICv2": "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V2/nCoV-2019.bed",
        "SARS-CoV-2-ARTICv3": "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V3/nCoV-2019.bed",
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
        write_artic_style_bed(StringIO(response.text), bed_output_filename)
        description = name[:-2] + " " + name[-2:] + " primer set"
        data.append(dict(value=name, path=bed_output_filename, description=description))
    return data


def install_primer_file(
    output_directory, input_filename, primer_name, primer_description
):
    name = re.sub(r"\W", "", str(primer_name).replace(" ", "_"))
    output_filename = os.path.join(output_directory, name + ".bed")
    with open(input_filename) as input_file:
        write_artic_style_bed(input_file, output_filename)
    data = [dict(value=name, description=primer_description, path=output_filename)]
    return data


class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(","))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch ARTIC SARS-CoV-2 primer files for Galaxy/IRIDA use"
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
    data_manager_dict["data_tables"][DATA_TABLE_NAME] = data_manager_dict[
        "data_tables"
    ].get(DATA_TABLE_NAME, [])

    if args.artic_primers:
        data = fetch_artic_primers(output_directory, args.artic_primers)
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
