#!/usr/bin/env python

from __future__ import print_function, division

import argparse
import json
import os
import os.path
import sys

import requests

DATA_TABLE_NAME = "artic_sars_cov2_primers"


def fetch_artic_primers(output_filename, output_directory, primers):
    primer_sets = {
        "ARTICv1": "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V1/nCoV-2019.bed",
        "ARTICv2": "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V2/nCoV-2019.bed",
        "ARTICv3": "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V3/nCoV-2019.bed",
    }

    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)
    data_manager_dict = {}
    data_manager_dict["data_tables"] = json.load(open(output_filename)).get(
        "data_tables", {}
    )
    data_manager_dict["data_tables"] = data_manager_dict.get("data_tables", {})
    data_manager_dict["data_tables"][DATA_TABLE_NAME] = data_manager_dict[
        "data_tables"
    ].get(DATA_TABLE_NAME, [])

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
        open(bed_output_filename, "w").write(response.text)
        description = name[:-2] + " " + name[-2:] + " primer set"
        data.append(dict(value=name, path=bed_output_filename, description=description))
    print(data)
    data_manager_dict["data_tables"][DATA_TABLE_NAME].extend(data)
    print(data_manager_dict)
    json.dump(data_manager_dict, open(output_filename, "w"))


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
    parser.add_argument(
        "--galaxy_datamanager_filename",
        help="Galaxy JSON format file describing data manager inputs",
    )
    parser.add_argument(
        "--primers",
        default="ARTICv1,ARTICv2,ARTICv3",
        action=SplitArgs,
        help="Comma separated list of primers to fetch",
    )
    args = parser.parse_args()

    config = json.load(open(args.galaxy_datamanager_filename))
    output_directory = config.get("output_data", [{}])[0].get("extra_files_path", None)
    if output_directory is None:
        output_directory = args.output_directory
    fetch_artic_primers(
        args.galaxy_datamanager_filename, output_directory, args.primers
    )
