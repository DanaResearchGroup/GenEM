from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import random
import logging
from molecules import Molecule


class MolecularSpace:
    """
    Class MolecularSpace
    Argument:
        initial_population(list): List of initial molecules
    """
    def __init__(self, initial_population):
        self.population = [Molecule(smiles) for smiles in initial_population]

    def calculate_distance(self, mol1, mol2):
        """
        Function to calculate the distance between two molecules
        Argument:
            self(Molecule): Molecule
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

    def mutate_molecule(self, molecule):
        """
        Function to mutate a molecule. It works by going through the list of mutations, randomly picking one. If the
        mutation results in a valid molecule, it returns the mutated molecule. Otherwise, it will randomly pick another
        mutation, and repeat the process. If the mutation doesn't produce a valid molecule after 3 times,
        it returns the original molecule.
        Argument:
            self(Molecule): Molecule
            molecule(Molecule): Molecule which is represented as fingerprint
        Return:
            mutated_molecule(Molecule): Mutated molecule
        """
        mutation_strategies = [
            self.smart_atom_substitution,
            self.isostere_replacement,
            self.functional_group_addition
        ]
        for _ in range(3):
            strategy = random.choice(mutation_strategies)
            mutated_molecule = strategy(molecule)
            if mutated_molecule and mutated_molecule.mol:
                return mutated_molecule
        return molecule

    def smart_atom_substitution(self, molecule):
        """
        Function to substitute atom with smart atom in order to increase chance of molecular validity.
        It defines allowed substitutions based on atom type and valence
        Argument:
            self(Molecule): Molecule
            molecule(Molecule): Molecule which is represented as fingerprint
        Return:
            mol(Molecule): Molecule which is represented as fingerprint
        """
        mol = Chem.RWMol(molecule.mol)
        if mol.GetNumAtoms() > 0:
            atom_idx = random.randint(0, mol.GetNumAtoms() - 1)
            atom = mol.GetAtomWithIdx(atom_idx)

            substitutions = {
                'C': ['N', 'O', 'S'],
                'N': ['C', 'O', 'S'],
                'O': ['N', 'S'],
                'S': ['O', 'N']
            }

            original_atom_symbol = atom.GetSymbol()
            if original_atom_symbol in substitutions:
                new_atom_symbol = random.choice(substitutions[original_atom_symbol])
                new_atom = Chem.Atom(new_atom_symbol)
                new_atom.SetFormalCharge(0)  # Reset formal charge
                mol.ReplaceAtom(atom_idx, new_atom)

        return self.sanitize_and_optimize_molecule(mol)

    from rdkit import Chem
    from rdkit.Chem import AllChem

    def find_connectable_atoms(self, mol):
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

    def combine_fragments(self, smiles):
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
                combined_connectable_atoms = self.find_connectable_atoms(combined_mol)
                frag_connectable_atoms = self.find_connectable_atoms(frag_conf)

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

    def isostere_replacement(self, molecule):
        """
        Function to replace similar functional groups with other functional groups.
        Argument:
        self(Molecule): Molecule
        molecule(Molecule): Molecule which is represented as fingerprint
        Return:
        mol(Molecule): Molecule which is represented as fingerprint
        """
        original_smiles = Chem.MolToSmiles(molecule.mol)
        mol = Chem.RWMol(molecule.mol)
        isosteres = {
            "C(=O)O": "C(=O)NH2",  # Carboxylic acid to amide
            "C(F)(F)F": "C#N",  # Trifluoromethyl to nitrile
            "c1ccccc1": "c1ccncc1",  # Benzene to pyridine
            "C1=CC2=C(C=C1)C(=O)N2C(=O)O": "C1=CC2=C(C=C1)C(=O)N2C(=O)O",
            # Nitro group replacement for energetic materials
            "C[N+](=O)[O-]": "C1=CC=C(C(=O)O)C=C1",  # Nitro to carboxyl group
            # Nitro group replacement to fused ring system
            "[O-][N+](=O)C": "[O-][N+](=O)C1=CC=CC=C1",  # Nitro group to phenyl ring
            "[N+](=O)[O-]": "[N+](=O)[O-]C1=CC=C(C(=O)O)C=C1"  # Nitro group to substituted phenyl ring
        }

        try:
            for original, replacement in isosteres.items():
                patt = Chem.MolFromSmarts(original)
                repl = Chem.MolFromSmiles(replacement)

                if patt and repl and mol.HasSubstructMatch(patt):
                    mol = AllChem.ReplaceSubstructs(mol, patt, repl, replaceAll=True)[0]
                    break
        except Exception:
            return None  # Return None if the replacement fails

        smiles = Chem.MolToSmiles(mol)

        # Check if the molecule was changed
        if smiles == original_smiles:
            return None

        # Combine fragments if necessary
        if '.' in smiles:
            try:
                combined_mol = self.combine_fragments(smiles)
                if combined_mol:
                    mol = combined_mol
            except Exception:
                return None

        return self.sanitize_and_optimize_molecule(mol)

    def functional_group_addition(self, molecule):
        """
        Function to add functional groups to molecule.
        Argument:
            self(Molecule): Molecule
            molecule(Molecule): Molecule which is represented as fingerprint
        Return:
            mol(Molecule): Molecule which is represented as fingerprint
        """
        mol = Chem.RWMol(molecule.mol)
        functional_groups = ["[OH]", "[NH2]", "[C](=O)[OH]", "[CH3]"]
        if mol.GetNumAtoms() > 0:
            atom_idx = random.randint(0, mol.GetNumAtoms() - 1)
            fg = Chem.MolFromSmarts(random.choice(functional_groups))
            if fg:
                edmol = Chem.EditableMol(mol)
                edmol.ReplaceAtom(atom_idx, fg.GetAtomWithIdx(0))
                mol = edmol.GetMol()

        return self.sanitize_and_optimize_molecule(mol)


    def sanitize_and_optimize_molecule(self, mol):
        """
        Function to sanitize the molecule after mutation to improve validity. In addition, logging has been designed
        to track invalid molecules and what mutations are causing them.
        Argument:
            self(Molecule): Molecule which is represented as fingerprint
            mol(Molecule): Molecule which is represented as fingerprint after mutation
        Return:
            Molecule(new_smiles): Molecule which is represented as fingerprint after undergoing sanitization.
        """
        try:
            Chem.SanitizeMol(mol)
            mol = Chem.AddHs(mol)
            new_smiles = Chem.MolToSmiles(mol)
            return Molecule(new_smiles)
        except Exception:
            return None
