from __future__ import annotations

import random

from rdkit import Chem
from rdkit.Chem import AllChem

from src.helpers.constants import ISOSTERES
from src.helpers.molecule_helpers import combine_fragments


class Molecule:
    """
    Class Molecule
    Args:
        smiles (str): SMILES string
        initial_population(list): List of initial molecules
    """

    def __init__(self, smiles):
        self.smiles = smiles
        self.properties = {}
        self.best_solution = None
        self.best_fitness = -float("inf")

        try:
            self.rdkit_mol = Chem.MolFromSmiles(smiles)
            if self.rdkit_mol is None:
                raise ValueError(f"Invalid SMILES string: {smiles}")
        except Exception as e:
            print(f"Error: {e}")
            self.rdkit_mol = None  # Set to None if SMILES parsing fails

    @staticmethod
    def smart_atom_substitution(mol):
        """
        Function to substitute atom with smart atom in order to increase chance of molecular validity.
        It defines allowed substitutions based on atom type and valence
        Argument:
            self(Molecule): Molecule
            molecule(Molecule): Molecule which is represented as fingerprint
        Return:
            mol(Molecule): Molecule which is represented as fingerprint
        """
        mol = Chem.RWMol(mol.rdkit_mol)
        if mol.GetNumAtoms() > 0:
            atom_idx = random.randint(0, mol.GetNumAtoms() - 1)
            atom = mol.GetAtomWithIdx(atom_idx)

            SUBSTITUTIONS = {
                "C": ["N", "O", "S"],
                "N": ["C", "O", "S"],
                "O": ["N", "S"],
                "S": ["O", "N"],
            }

            original_atom_symbol = atom.GetSymbol()
            if original_atom_symbol in SUBSTITUTIONS:
                new_atom_symbol = random.choice(SUBSTITUTIONS[original_atom_symbol])
                new_atom = Chem.Atom(new_atom_symbol)
                new_atom.SetFormalCharge(0)  # Reset formal charge
                mol.ReplaceAtom(atom_idx, new_atom)

        return Molecule.sanitize_and_optimize_molecule(mol)

    @staticmethod
    def isostere_replacement(origin_mol: Molecule) -> Molecule | None:
        """
        Function to replace similar functional groups with other functional groups.
        Argument:
        self(Molecule): Molecule
        molecule(Molecule): Molecule which is represented as fingerprint
        Return:
        mol(Molecule): Molecule which is represented as fingerprint
        """

        mutated_mol = Chem.RWMol(origin_mol.rdkit_mol)

        try:
            for pattern, replacement in ISOSTERES.items():
                pattern_mol = Chem.MolFromSmiles(pattern)
                replacement_mol = Chem.MolFromSmiles(replacement)

                if (
                    pattern_mol
                    and replacement_mol
                    and origin_mol.rdkit_mol.HasSubstructMatch(pattern_mol)
                ):
                    all_possible_replacements = AllChem.ReplaceSubstructs(
                        mutated_mol, pattern_mol, replacement_mol
                    )

                    # randomly choose only one of the replacements
                    mutated_mol = all_possible_replacements[
                        random.randint(0, len(all_possible_replacements) - 1)
                    ]
                    break
        except Exception:
            # TODO: add log or handler for this exception
            return None  # Return None if the replacement fails

        new_mol_smiles = Chem.MolToSmiles(mutated_mol)

        # Check if the molecule was changed
        if new_mol_smiles == origin_mol.properties.get("smiles"):
            return None

        # Combine fragments if necessary
        if "." in new_mol_smiles:
            try:
                combined_mol = combine_fragments(new_mol_smiles)
                if combined_mol:
                    mutated_mol = combined_mol
            except Exception:
                return None

        return Molecule.sanitize_and_optimize_molecule(mutated_mol)

    @staticmethod
    def functional_group_addition(mol):
        """
        Function to add functional groups to molecule.
        Argument:
            self(Molecule): Molecule
            molecule(Molecule): Molecule which is represented as fingerprint
        Return:
            mol(Molecule): Molecule which is represented as fingerprint
        """
        mol = Chem.RWMol(mol.rdkit_mol)
        FUNCTIONAL_GROUPS = ["[OH]", "[NH2]", "[C](=O)[OH]", "[CH3]"]
        if mol.GetNumAtoms() > 0:
            atom_idx = random.randint(0, mol.GetNumAtoms() - 1)
            fg = Chem.MolFromSmarts(random.choice(FUNCTIONAL_GROUPS))
            if fg:
                edmol = Chem.EditableMol(mol)
                edmol.ReplaceAtom(atom_idx, fg.GetAtomWithIdx(0))
                mol = edmol.GetMol()

        return Molecule.sanitize_and_optimize_molecule(mol)

    @staticmethod
    def sanitize_and_optimize_molecule(mol: Molecule) -> Molecule | None:
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
