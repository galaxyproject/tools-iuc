#!/usr/bin/env python

import argparse
import json
import os
import sys
from pathlib import Path

import requests

# URL for downloading OMAmer datasets
OMAMER_DATASETS_URL = "https://omabrowser.org/All/{dataset}"


# List of OMAmer data sets
OMAMER_DATASETS = [
    "Primates-v2.0.0.h5",
    "Viridiplantae-v2.0.0.h50",
    "Metazoa-v2.0.0.h5",
    "LUCA-v0.2.5.h5"
]


def download_file(url, dest):
    try:
        # Download file from URL to local location
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(dest, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        print(f"Downloaded: {url} to {dest}")
    except requests.exceptions.RequestException as e:
        # Handles download errors
        print(f"Error downloading {url}: {e}")
        sys.exit(1)


def main(args):
    # Create output directory if none exists
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)


    # Download each specified OMAmer data set
    for dataset in OMAMER_DATASETS:
        url = OMAMER_DATASETS_URL.format(dataset=dataset)
        base_name = os.path.splitext(dataset)[0]
        destination_path = os.path.join(args.output_dir, base_name)
        download_file(url, destination_path)


    # Utiliser le nom du fichier sans extension comme identifiant unique
    data_manager_entry = {
        "value": os.path.splitext(os.path.basename(base_name))[0],
        "name": os.path.splitext(os.path.basename(base_name))[0],
        "path": str(Path(args.output_dir)),
    }

    # Creates a JSON dictionary representing the Data Manager configuration
    data_manager_json = {"data_tables": {"omamer_data": [data_manager_entry]}}


    # Writes this JSON dictionary to the specified output file
    with open(args.json, "w") as fh:
        json.dump(data_manager_json, fh, indent=2, sort_keys=True)


if __name__ == "__main__":
    # Set up argparse to specify expected command line arguments
    parser = argparse.ArgumentParser(description='Download data for OMAmer')
    parser.add_argument('--output-dir', dest='output_dir', required=True, help='Output directory for saving databases')
    parser.add_argument('--name', default='default_name', help='Data table entry unique ID')
    parser.add_argument('--json', help='Path to JSON file')
    args = parser.parse_args()


    # Call the main function with the analyzed arguments
    main(args)