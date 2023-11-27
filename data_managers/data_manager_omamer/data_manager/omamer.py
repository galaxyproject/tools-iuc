#!/usr/bin/env python

import argparse
import json
import os
import subprocess
from pathlib import Path
from datetime import datetime
import sys

import requests

OMAMER_DATASETS_URL = "https://omabrowser.org/All/{dataset}"
OMAMER_DATASETS = [
    "Homininae.h5",
    "Primates-v2.0.0.h5",
    "Saccharomyceta.h5",
    "Viridiplantae-v0.2.5.h50",
    "Metazoa-v0.2.5.h5",
    "LUCA-v0.2.5.h5"
]

def download_file(url, dest):
    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(dest, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        print(f"Downloaded: {url} to {dest}")
    except requests.exceptions.RequestException as e:
        print(f"Error downloading {url}: {e}")
        sys.exit(1)

def main(args):
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    for dataset in OMAMER_DATASETS:
        url = OMAMER_DATASETS_URL.format(dataset=dataset)
        base_name = os.path.splitext(dataset)[0]
        destination_path = os.path.join(args.output_dir, base_name)
        download_file(url, destination_path)

    data_manager_entry = {
        "value": args.name.lower(),
        "name": args.name,
        "path": str(Path(args.output_dir)),
    }
    data_manager_json = {"data_tables": {"omamer_data": [data_manager_entry]}}

    with open(args.json, "w") as fh:
        json.dump(data_manager_json, fh, indent=2, sort_keys=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download data for OMAmer')
    parser.add_argument('--output-dir', dest='output_dir', required=True, help='Output directory for saving databases')
    parser.add_argument('--name', default=str(datetime.date.today()), help='Data table entry unique ID')
    parser.add_argument('--json', help='Path to JSON file')
    args = parser.parse_args()

    main(args)