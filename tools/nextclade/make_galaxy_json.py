#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate galaxy.json to describe nextclade column names"
    )
    parser.add_argument(
        "--output_name",
        default="report_tsv",
        help="Name of the Galaxy output to add metadata for",
    )
    parser.add_argument(
        "nextclade_output_file",
        type=argparse.FileType("r"),
        help="Nextclade output in TSV format",
    )
    args = parser.parse_args()

    headers = next(args.nextclade_output_file).strip().split("\t")
    galaxy_json = {args.output_name: {"metadata": {"column_names": headers}}}
    print(json.dumps(galaxy_json, indent=2))
