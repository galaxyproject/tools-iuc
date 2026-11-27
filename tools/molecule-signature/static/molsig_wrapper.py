#!/usr/bin/env python
import json
import subprocess
import sys

from molsig.enumerate_utils import get_first_stereoisomer
from molsig.SignatureAlphabet import SignatureAlphabet, load_alphabet, merge_alphabets
from rdkit import Chem

RADIUS = 2
NBITS = 2048
USE_STEREO = True


def main():
    inputs_json, alphabet_path = sys.argv[1], sys.argv[2]
    with open(inputs_json) as fh:
        params = json.load(fh)

    smiles_0 = params["input_text"].strip()
    include_input = params["incl_input_alf"]

    if Chem.MolFromSmiles(smiles_0) is None:
        sys.exit(f"SMILES is not validated by rdkit: {smiles_0}")

    if include_input:
        alphabet = load_alphabet(alphabet_path)
        smiles = get_first_stereoisomer(smiles_0)
        alphabet_mol = SignatureAlphabet(radius=RADIUS, nBits=NBITS, use_stereo=USE_STEREO)
        alphabet_mol.fill([smiles], verbose=False)
        alphabet = merge_alphabets(alphabet, alphabet_mol)
        alphabet_file = "updated_alphabet.npz"
        alphabet.save(alphabet_file)
    else:
        alphabet_file = alphabet_path

    result = subprocess.run([
        "molsig", "enumerate",
        f"--smiles={smiles_0}",
        "--alphabet", alphabet_file,
        "--output", "analogues_smiles.tsv",
    ])
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()