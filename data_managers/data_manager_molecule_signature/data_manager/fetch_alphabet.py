#!/usr/bin/env python
import argparse
import hashlib
import json
import logging
import os
import shutil
import sys
import urllib.request

ZENODO = "https://zenodo.org/records/22691018/files"

ALPHABETS = {
    "chembl": {
        "name": "ChEMBL alphabet",
        "file": "chembl_alphabet.npz",
        "md5": "9c1f00da7f1bc300c2b7682f9bc06ed0",
    },
    "drugbank": {
        "name": "DrugBank alphabet",
        "file": "drugbank_alphabet.npz",
        "md5": "62fedf7cceabb6553dec0efd58f0f763",
    },
    "emolecules": {
        "name": "eMolecules alphabet",
        "file": "emolecules_alphabet.npz",
        "md5": "36fccb3f550e94713a7764e92a086481",
    },
    "metanetx": {
        "name": "MetaNetX alphabet",
        "file": "metanetx_alphabet.npz",
        "md5": "0eac175891133b29d9172fc14c487ae2",
    },
}


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stderr)],
    )
    parser = argparse.ArgumentParser()
    parser.add_argument("--alphabet", required=True, choices=ALPHABETS)
    parser.add_argument("--version", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    # The pre-existing "out" file holds the JSON mapping of the job
    # parameters; the directory where the alphabet file is staged for
    # Galaxy to move into the managed data directory is read from the
    # first output's extra_files_path.
    with open(args.out) as fh:
        params = json.load(fh)
    target_dir = params["output_data"][0]["extra_files_path"]
    os.makedirs(target_dir, exist_ok=True)

    info = ALPHABETS[args.alphabet]
    url = f"{ZENODO}/{info['file']}?download=1"
    dest = os.path.join(target_dir, info["file"])
    tmp = dest + ".part"

    logging.info("Downloading %s", url)
    with urllib.request.urlopen(url) as resp, open(tmp, "wb") as fh:
        shutil.copyfileobj(resp, fh)
    logging.info("Downloaded %d bytes", os.path.getsize(tmp))

    observed = md5sum(tmp)
    if observed != info["md5"]:
        os.remove(tmp)
        logging.error(
            "MD5 mismatch for %s: expected %s, got %s",
            info["file"],
            info["md5"],
            observed,
        )
        sys.exit(1)
    os.replace(tmp, dest)
    logging.info("MD5 verified for %s", info["file"])

    entry = {
        "value": args.alphabet,
        "name": info["name"],
        "version": args.version,
        # Relative to extra_files_path; data_manager_conf.xml moves the file
        # from there into the managed data directory.
        "path": info["file"],
    }
    with open(args.out, "w") as fh:
        json.dump({"data_tables": {"molsig_alphabets": [entry]}}, fh, sort_keys=True)


if __name__ == "__main__":
    main()
