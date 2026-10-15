#!/usr/bin/env python

import argparse
import hashlib
import json
import os
import sys
import urllib.request
import zipfile

URL = "https://zenodo.org/records/10041121/files/data.zip?download=1" 


def sha256_of(path, chunk_size=1024 * 1024):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, help="Directory to create and fill with the data")
    parser.add_argument("--output-file", required=True, help="Path to JSON file")
    parser.add_argument("--test", default="false", help="'true' to skip download and create an empty data dir")
    args = parser.parse_args()

    output_dir = os.path.abspath(args.output_dir)
    output_file = os.path.abspath(args.output_file)
    test = args.test.strip().lower() in ("true", "1", "yes")

    os.mkdir(output_dir)
    os.chdir(output_dir)

    checksum = ""
    if test:
        os.mkdir("data")
    else:
        urllib.request.urlretrieve(URL, "data.zip")
        checksum = sha256_of("data.zip")
        with zipfile.ZipFile("data.zip") as z:
            z.extractall(".")

    if not os.path.isdir("data"):
        sys.stderr.write("data missing\n")
        sys.exit(1)

    dm = {
        "data_tables": {
            "enzbert": [
                {
                    "value": "EnzBert_SwissProt_2021_04",
                    "name": "Enzbert_SwissProt_2021_04",
                    "path": output_dir,
                    "checksum": checksum,
                },
                {
                    "value": "EnzBert_EC40",
                    "name": "Enzbert__EC40",
                    "path": output_dir,
                    "checksum": checksum,
                },
                {
                    "value": "EnzBert_ECPred40",
                    "name": "Enzbert_ECPred40" ,
                    "path": output_dir,
                    "checksum": checksum,
                },
            ]
        }
    }

    with open(output_file, "w") as f:
        json.dump(dm, f, indent=2)


if __name__ == "__main__":
    main()