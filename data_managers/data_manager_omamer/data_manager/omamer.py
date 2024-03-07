#!/usr/bin/env python

import argparse
import json
import os
import sys
from pathlib import Path

import requests

# URL for downloading OMAmer datasets
OMAMER_URL = "https://omabrowser.org/All/{dataset}"

# List of OMAmer data sets with versions
OMAMER = {
    "Primates": "Primates-v2.0.0.h5",
    "Viridiplantae": "Viridiplantae-v2.0.0.h5",
    "Metazoa": "Metazoa-v2.0.0.h5",
    "LUCA": "LUCA-v0.2.5.h5",
}

DEFAULT_OUTPUT_DIR = "database_omamer"


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
    # Set output directory to default
    output_dir = DEFAULT_OUTPUT_DIR

    # Create output directory if none exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Check if the selected dataset exists
    if args.name not in OMAMER:
        print(f"Error: Selected dataset '{args.name}' not found.")
        sys.exit(1)

    # Download the selected OMAmer dataset
    dataset = OMAMER[args.name]
    url = OMAMER_URL.format(dataset=dataset)
    base_name = os.path.splitext(dataset)[0]
    destination_path = os.path.join(output_dir, base_name)
    download_file(url, destination_path)

    data_manager_entry = {
        "value": os.path.splitext(os.path.basename(base_name))[0],
        "name": os.path.splitext(os.path.basename(base_name))[0],
        "version": "2.0.2",
        "path": str(Path(output_dir)),
    }

    # Creates a JSON dictionary representing the Data Manager configuration
    data_manager_json = {"data_tables": {"omamer_data": [data_manager_entry]}}

    # Writes this JSON dictionary to the specified output file
    with open(args.json, "w") as fh:
        json.dump(data_manager_json, fh, indent=2, sort_keys=True)


if __name__ == "__main__":
    # Set up argparse to specify expected command line arguments
    parser = argparse.ArgumentParser(description='Download data for OMAmer')
    parser.add_argument('--name', default='Primates', choices=OMAMER.keys(), help='Select dataset to download')
    parser.add_argument('--json', help='Path to JSON file')

    args = parser.parse_args()

    main(args)
