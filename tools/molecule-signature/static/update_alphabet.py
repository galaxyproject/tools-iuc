import sys
from rdkit import Chem
from molsig.enumerate_utils import get_first_stereoisomer
from molsig.SignatureAlphabet import (SignatureAlphabet,merge_alphabets,load_alphabet)

RADIUS = 2
NBITS = 2048
USE_STEREO = True

alphabet_path = sys.argv[1]
output_path   = sys.argv[2]
smiles_0        = sys.argv[3]

Alphabet = load_alphabet(alphabet_path)

mol = Chem.MolFromSmiles(smiles_0)
smiles = get_first_stereoisomer(smiles_0)

Alphabet_mol = SignatureAlphabet(radius=RADIUS, nBits=NBITS, use_stereo=USE_STEREO)
Alphabet_mol.fill([smiles], verbose=False)
Alphabet = merge_alphabets(Alphabet, Alphabet_mol)

Alphabet.save(output_path)