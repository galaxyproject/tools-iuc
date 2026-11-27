#!/usr/bin/env python
import argparse
import hashlib
import json
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
    parser = argparse.ArgumentParser()
    parser.add_argument("--alphabet", required=True, choices=ALPHABETS)
    parser.add_argument("--version", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    # Le fichier "out" contient les paramètres du job ; on y lit le dossier cible
    with open(args.out) as fh:
        params = json.load(fh)
    target_dir = params["output_data"][0]["extra_files_path"]
    os.makedirs(target_dir, exist_ok=True)

    info = ALPHABETS[args.alphabet]
    url = f"{ZENODO}/{info['file']}?download=1"
    dest = os.path.join(target_dir, info["file"])
    tmp = dest + ".part"

    with urllib.request.urlopen(url) as resp, open(tmp, "wb") as fh:
        shutil.copyfileobj(resp, fh)

    observed = md5sum(tmp)
    if observed != info["md5"]:
        os.remove(tmp)
        sys.exit(
            f"MD5 mismatch for {info['file']}: expected {info['md5']}, got {observed}"
        )
    os.replace(tmp, dest)

    entry = {
        "value": args.alphabet,
        "name": info["name"],
        "version": args.version,
        "path": dest,
    }
    with open(args.out, "w") as fh:
        json.dump({"data_tables": {"molsig_alphabets": [entry]}}, fh, sort_keys=True)


if __name__ == "__main__":
    main()
