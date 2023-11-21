#!/usr/bin/env python

import argparse
import json
import operator
import os
import subprocess
import sys
import tarfile
from datetime import datetime

import requests

OMAMER_DATASETS_URL = "https://omabrowser.org/All/{dataset}"
OMAMER_DATASETS = [
    "Homininae.h5",
    "Primates-v2.0.0.h5",
    "Saccharomyceta.h5",
    "Viridiplantae-v0.2.5.h50",
    "Viridiplantae-v2.0.0.h5",
    "Metazoa-v0.2.5.h5",
    "Metazoa-v2.0.0.h5",
    "LUCA.h5",
    "LUCA-v2.0.0.h5",
    "LUCA-v0.2.5.h5"
]

def download_file(url, dest):
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(dest, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', dest='output_dir', required=True, help='Output directory for saving databases')
    args = parser.parse_args()

    for dataset in OMAMER_DATASETS:
        url = OMAMER_DATASETS_URL.format(dataset=dataset)
        base_name = os.path.splitext(dataset)[0]
        destination_path = os.path.join(args.output_dir, base_name)
        download_file(url, destination_path)


