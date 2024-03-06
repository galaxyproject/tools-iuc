#!/usr/bin/env python

import argparse
import json
import os
import subprocess


def main():
    opts = parse_args()

    output_dict = {
        "data_tables": {
            "ncbi_fcs_gx_databases": sync_files(opts),
            "ncbi_fcs_gx_divisions": get_divisions(opts),
        }
    }

    with open(opts.output_file, "w") as f:
        print(json.dumps(output_dict, sort_keys=True, indent=2), file=f)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", required=True)
    parser.add_argument("--source_manifest", required=True)
    parser.add_argument("--output_file", required=True)
    parser.add_argument("--output_dir", required=True)

    return parser.parse_args()


def sync_files(opts):
    os.makedirs(opts.output_dir, exist_ok=True)

    args = [
        "sync_files.py",
        "--mft",
        opts.source_manifest,
        "--dir",
        opts.output_dir,
        "get",
    ]

    try:
        subprocess.run(args, capture_output=True, check=True)
    except subprocess.CalledProcessError:
        raise

    entries_dict = {
        "add": [
            {
                "value": opts.tag,
                "source_manifest": opts.source_manifest,
                "name": opts.output_dir,
            }
        ]
    }

    return entries_dict


def get_divisions(opts):
    # descriptions for the top-level gx divisions
    top_level_description = {
        "anml": "Animals (Metazoa)",
        "arch": "Archaea",
        "fung": "Fungi",
        "plnt": "Plants (Viridiplantae)",
        "prok": "Bacteria",
        "prst": "Protists (other Eukaryota)",
        "synt": "Synthetic",
        "virs": "Virus",
    }

    # get the pathname for the taxa file
    manifest_filename = os.path.basename(opts.source_manifest)
    assert manifest_filename.lower().endswith(
        ".manifest"
    ), 'source_manifest does not end with ".manifest"'
    manifest_tag = manifest_filename[:-9]
    taxa_pathname = os.path.join(opts.output_dir, f"{manifest_tag}.taxa.tsv")

    gx_divisions = set()
    with open(taxa_pathname) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            tax_id, species, common_name, blast_div, div = line.split("\t", 4)
            gx_divisions.add(div)

    elements = []
    for division in gx_divisions:
        top, bottom = division.split(":", 1)
        description = f"{top_level_description[top]} - {bottom}"
        elements.append((description, division))

    # add an element to support unknown/unclassified samples
    elements.append(("Unknown / Unclassified", "unkn:unknown"))

    entries_dict = {"add": []}

    for name, gx_div in sorted(elements):
        entries_dict["add"].append({"value": gx_div, "tag": opts.tag, "name": name})

    return entries_dict


if __name__ == "__main__":
    main()
