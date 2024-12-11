from rdkit import RDLogger

from src.helpers.constants import MOLECULES_SMILES
from src.molecule import Molecule

# Set log level to suppress warnings
RDLogger.logger().setLevel(RDLogger.CRITICAL)

INITIAL_SMILES = [
    MOLECULES_SMILES["Aspirin"],
    MOLECULES_SMILES["Testosterone"],
    MOLECULES_SMILES["Nicotine"],
    MOLECULES_SMILES["TNT"],
    MOLECULES_SMILES["RDX"],
    MOLECULES_SMILES["HMX"],
    MOLECULES_SMILES["PETN"],
    MOLECULES_SMILES["Tetryl"],
    # TODO: find out names (if any) and add to MOLECULES_SMILES
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "CC(=O)NC1=CC=C(O)C=C1",
    "CN(C)CCN1C2=CC=CC=C2SC3=CC=CC=C13",
    "CN(C)C(=O)C1=CC2=C(N1)C=CC3=C2C=CC4=C3C=CN4C",
]
count = 0
hist = [[0, set()] for _ in INITIAL_SMILES]
for _ in range(20):
    for j in range(len(INITIAL_SMILES)):
        mutated_mol = Molecule.functional_group_addition(Molecule(INITIAL_SMILES[j]))
        if mutated_mol:
            hist[j][0] += 1
            hist[j][1].add(mutated_mol.smiles)
            count += 1

print(count, "out of", len(INITIAL_SMILES) * 20)
molecule_keys_list = list(MOLECULES_SMILES.keys())
for h in range(len(hist)):
    print(molecule_keys_list[h],":", INITIAL_SMILES[h], ",",hist[h][0])
    for j in hist[h][1]:
        print(j)



