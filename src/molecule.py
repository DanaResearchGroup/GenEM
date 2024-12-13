from __future__ import annotations

import random

from rdkit import Chem
from rdkit.Chem import AllChem

from src.helpers.constants import BACKBONE_LIST, FUNCTIONAL_GROUPS, SUBSTITUTIONS
from src.helpers.molecule_helpers import combine_fragments


class Molecule:
    """
    Class Molecule
    Args:
        smiles (str): SMILES string
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
    def smart_atom_substitution(origin_mol: Molecule) -> Molecule | None:
        """
        Function to substitute atom with another atom.
        Argument:
            mol(Molecule): Molecule class instance
        Return:
            mol(Molecule): Molecule class instance
        """
        mutated_mol = Chem.RWMol(origin_mol.rdkit_mol)
        if mutated_mol.GetNumAtoms() > 0:
            atom_idx = random.randint(0, mutated_mol.GetNumAtoms() - 1)
            atom = mutated_mol.GetAtomWithIdx(atom_idx)
            original_atom_symbol = atom.GetSymbol()
            if original_atom_symbol in SUBSTITUTIONS:
                new_atom_symbol = random.choice(SUBSTITUTIONS)
                new_atom = Chem.Atom(new_atom_symbol)
                new_atom.SetFormalCharge(0)  # Reset formal charge
                mutated_mol.ReplaceAtom(atom_idx, new_atom)

        return Molecule.sanitize_and_optimize_molecule(mutated_mol)

    @staticmethod
    def backbone_replacement(origin_mol: Molecule) -> Molecule | None:
        """
        Function to replace common backbones  with other backbones often found in EMs.
        Argument:
        origin_mol(Molecule): Molecule from class Molecule
        Return:
        mol(Molecule): Molecule from class Molecule
        """

        mutated_mol = Chem.RWMol(origin_mol.rdkit_mol)

        possible_matching_pattern = None

        # loop through all patterns to find optional replacements
        for pattern in random.sample(BACKBONE_LIST, len(BACKBONE_LIST)):
            pattern_mol = Chem.MolFromSmiles(pattern)
            if pattern_mol and origin_mol.rdkit_mol.HasSubstructMatch(
                    pattern_mol
            ):  # TODO: catch not valid 'pattern_mol' and log them
                possible_matching_pattern = pattern_mol
                # once found an optional patter, break the loop
                break

        # find all the possible replacements for a given pattern

        for optional_replacement in random.sample(BACKBONE_LIST, len(BACKBONE_LIST)):
            replacement_mol = Chem.MolFromSmiles(optional_replacement)
            if (
                    possible_matching_pattern
                    and replacement_mol
                    and optional_replacement != Chem.MolToSmiles(possible_matching_pattern)
            ):  # TODO: catch not valid 'replacement_mol' and log them
                possible_replacements = AllChem.ReplaceSubstructs(
                    mutated_mol, possible_matching_pattern, replacement_mol
                )

                if possible_replacements and len(possible_replacements):
                    mutated_mol = random.choice(possible_replacements)
                break

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
                # TODO: add log here
                return None

        return Molecule.sanitize_and_optimize_molecule(mutated_mol)

    @staticmethod
    def functional_group_addition(mol):
        """
        Function to add functional groups to molecule.
        Argument:
            mol(Molecule): molecule from class Molecule
        Return:
            mol(Molecule): molecule from class Molecule
        """
        mol = Chem.RWMol(mol.rdkit_mol)
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
        Function to sanitize the molecule after mutation to improve validity.
        Argument:
            mol(Molecule): Molecule from class Molecule
        Return:
            Molecule(new_smiles): Molecule from class Molecule after undergoing sanitization.
        """
        try:
            Chem.SanitizeMol(mol)
            mol = Chem.AddHs(mol)
            new_smiles = Chem.MolToSmiles(mol)
            return Molecule(new_smiles)
        except Exception:
            # TODO: log error here
            return None
