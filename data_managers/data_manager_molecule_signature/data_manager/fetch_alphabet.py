#!/usr/bin/env python
"""Download a molecule-signature alphabet file into a target directory.

The data manager JSON that describes the new data table entry is generated
by a configfile in the tool wrapper, so this script only downloads data.
"""
import argparse
import hashlib
import logging
import os
import shutil
import sys
import urllib.request

ZENODO = "https://zenodo.org/records/22691018/files"

# Upstream file names and checksums. The file names are repeated in the
# "dmjson" configfile of molsig_alphabet_fetcher.xml, which builds the data
# table row; keep the two in sync.
ALPHABETS = {
    "chembl": {
        "file": "chembl_alphabet.npz",
        "md5": "9c1f00da7f1bc300c2b7682f9bc06ed0",
    },
    "drugbank": {
        "file": "drugbank_alphabet.npz",
        "md5": "62fedf7cceabb6553dec0efd58f0f763",
    },
    "emolecules": {
        "file": "emolecules_alphabet.npz",
        "md5": "36fccb3f550e94713a7764e92a086481",
    },
    "metanetx": {
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--alphabet", required=True, choices=sorted(ALPHABETS), help="Alphabet to download"
    )
    parser.add_argument(
        "--out", required=True, help="Directory to download the alphabet file into"
    )
    args = parser.parse_args()

    info = ALPHABETS[args.alphabet]
    url = f"{ZENODO}/{info['file']}?download=1"
    os.makedirs(args.out, exist_ok=True)
    dest = os.path.join(args.out, info["file"])
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


if __name__ == "__main__":
    main()
