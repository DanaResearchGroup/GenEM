from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def calculate_distance(mol1, mol2):
    """
    Function to calculate the distance between two molecules
    Argument:
        mol1(Molecule): Molecule which is represented as fingerprint
        mol2(Molecule): Molecule which is represented as fingerprint
    Return:
        distance(float): Distance between two molecules
    """
    if mol1.fingerprint is None:
        mol1.calculate_fingerprint()
    if mol2.fingerprint is None:
        mol2.calculate_fingerprint()
    return 1 - DataStructs.TanimotoSimilarity(mol1.fingerprint, mol2.fingerprint)

def get_nearest_neighbors(self, mol, k=5):
    """
    Function to calculate the nearest neighbors of a given molecule
    Argument:
        self(Molecule): Molecule
        mol(Molecule): Molecule which is represented as fingerprint
        k(int): Number of nearest neighbors
    Return:
        neighbors(list): List of nearest neighbors
    """
    distances = [(m, self.calculate_distance(mol, m)) for m in self.population if m != mol]
    return sorted(distances, key=lambda x: x[1])[:k]

def find_connectable_atoms(mol):
    """
    Find atoms with unfilled valences (often terminal atoms) in a molecule.
    This helps to identify which atoms can be used to connect fragments.
    """
    connectable_atoms = []
    for atom in mol.GetAtoms():
        # Check if atom has unfilled valence (potential for bonding)
        if atom.GetNumImplicitHs() > 0:
            connectable_atoms.append(atom.GetIdx())
    return connectable_atoms

def combine_fragments(smiles):
        """
        Combine molecular fragments by connecting them with single bonds.
        Argument:
        - smiles (str): SMILES representation of fragmented molecules.
        Return:
        - combined_mol: Combined molecule as RDKit Mol object, or None if failed.
        """
        fragments = smiles.split('.')
        if len(fragments) <= 1:
            return Chem.MolFromSmiles(smiles)  # No fragments to combine

        frag_mols = [Chem.MolFromSmiles(frag) for frag in fragments if frag]

        if not frag_mols or frag_mols[0] is None:
            return None  # Return None if no valid fragments

        try:
            # Start with the first fragment as the base
            combined_mol = Chem.RWMol(frag_mols[0])
        except Exception:
            return None  # Return None if initialization fails

        for frag in frag_mols[1:]:
            try:
                if frag is None:
                    continue

                frag_conf = Chem.RWMol(frag)
                frag_conf.UpdatePropertyCache()

                # Find connectable atoms in both the combined molecule and current fragment
                combined_connectable_atoms = find_connectable_atoms(combined_mol)
                frag_connectable_atoms = find_connectable_atoms(frag_conf)

                if not combined_connectable_atoms or not frag_connectable_atoms:
                    continue  # Skip if no suitable atoms found to connect

                # Connect the first found connectable atom from combined_mol with frag
                atom_idx_1 = combined_connectable_atoms[0]
                atom_idx_2 = frag_connectable_atoms[0]

                # Combine the molecules by adding the fragment
                combined_mol.InsertMol(frag_conf)
                combined_mol.AddBond(atom_idx_1, combined_mol.GetNumAtoms() - frag_conf.GetNumAtoms() + atom_idx_2,
                                     Chem.BondType.SINGLE)

            except Exception:
                continue  # Skip this fragment and move on

        try:
            final_mol = combined_mol.GetMol()
            Chem.SanitizeMol(final_mol)  # Ensure the final molecule is valid
            return final_mol
        except Chem.SanitizeException:
            return None
