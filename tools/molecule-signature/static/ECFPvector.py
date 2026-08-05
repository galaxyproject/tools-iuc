import sys

from molsig.enumerate_signature import enumerate_molecule_from_morgan
from molsig.SignatureAlphabet import load_alphabet

ecfp_text, alphabet_path, output_path = sys.argv[1], sys.argv[2], sys.argv[3]

for token in ("__ob__", "__cb__", "__oc__", "__cc__", "[", "]", "{", "}"):
    ecfp_text = ecfp_text.replace(token, "")
ecfp_vector = [int(x) for x in ecfp_text.replace(",", " ").split()]

alphabet = load_alphabet(alphabet_path)
Ssig, Smol, Nsig, thresholds_reached, computational_times = enumerate_molecule_from_morgan(
    ecfp_vector, alphabet
)

with open(output_path, "w") as fd:
    fd.write("SMILES\n")
    for smi in sorted(Smol):
        fd.write(smi + "\n")
