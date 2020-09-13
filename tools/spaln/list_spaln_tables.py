#!/usr/bin/env python3

import argparse
import shlex
import sys
from subprocess import run
from typing import TextIO


def find_common_ancestor_distance(
    taxon: str, other_taxon: str, taxonomy_db_path: str, only_canonical: bool
):
    canonical = "--only_canonical" if only_canonical else ""
    cmd_str = f"taxonomy_util -d {taxonomy_db_path} common_ancestor_distance {canonical} '{other_taxon}' '{taxon}'"
    cmd = shlex.split(cmd_str)
    proc = run(cmd, encoding="utf8", capture_output=True)
    return proc


def find_distances(gnm2tab_file: TextIO, taxon: str, taxonomy_db_path: str):
    cmd = ["taxonomy_util", "-d", taxonomy_db_path, "get_id", taxon]
    proc = run(cmd, capture_output=True, encoding="utf8")
    if "not found in" in proc.stderr:
        exit("Error: " + proc.stderr.strip())
    for line in gnm2tab_file:
        fields = line.split("\t")
        (species_code, settings, other_taxon) = map(lambda el: el.strip(), fields[:3])
        proc = find_common_ancestor_distance(taxon, other_taxon, taxonomy_db_path, True)
        ancestor_info = proc.stdout.rstrip()
        if proc.stderr != "":
            print("Warning:", other_taxon, proc.stderr.rstrip(), file=sys.stderr)
        else:
            proc = find_common_ancestor_distance(
                taxon, other_taxon, taxonomy_db_path, False
            )
            non_canonical_distance = proc.stdout.split("\t")[0]
            print(
                non_canonical_distance,
                ancestor_info,
                species_code,
                settings,
                other_taxon,
                sep="\t",
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find distance to common ancestor")
    parser.add_argument(
        "--taxonomy_db", required=True, help="NCBI Taxonomy database (SQLite format)"
    )
    parser.add_argument(
        "--gnm2tab_file",
        required=True,
        type=argparse.FileType(),
        help="gnm2tab file from spal",
    )
    parser.add_argument("taxon")
    args = parser.parse_args()

    find_distances(args.gnm2tab_file, args.taxon, args.taxonomy_db)
