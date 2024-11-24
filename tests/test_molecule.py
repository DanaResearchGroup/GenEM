#!/usr/bin/env python3
# encoding: utf-8

import unittest

from constants import MOLECULES_SMILES
from rdkit import Chem

from ... import (
    Molecule,
    functional_group_addition,
    isostere_replacement,
    smart_atom_substitution,
)


class TestMolecule(unittest.TestCase):

    def setUp(self):
        self.ethanol_smiles = MOLECULES_SMILES.get("Ethanol")  # Ethanol
        self.tnt_smiles = MOLECULES_SMILES.get("TNT")  # TNT
        self.invalid_smiles = "COCOCOCOCOC11##111OOO"  # invalid smile

        self.ethanol = Molecule(self.ethanol_smiles)
        self.tnt = Molecule(self.tnt_smiles)
        self.invalid_molecule = Molecule(self.invalid_smiles)

    def test_smiles_to_mol(self):
        expected_eth_smiles = Chem.CanonSmiles(
            self.ethanol_smiles
        )  # Canonicalize SMILES
        generated_eth_smiles = Chem.CanonSmiles(Chem.MolToSmiles(self.ethanol.mol))
        expected_tnt_smiles = Chem.CanonSmiles(self.tnt_smiles)  # Canonicalize SMILES
        generated_tnt_smiles = Chem.CanonSmiles(Chem.MolToSmiles(self.tnt.mol))
        self.assertEqual(generated_eth_smiles, expected_eth_smiles)
        self.assertEqual(generated_tnt_smiles, expected_tnt_smiles)

    def test_fingerprint(self):
        self.ethanol.calculate_fingerprint()
        self.assertIsNotNone(self.ethanol.fingerprint)
        self.assertEqual(len(self.ethanol.fingerprint), 1024)  # Default nBits=1024

    def _run_mutation_test(self, mutation_func):
        for mol in [self.ethanol, self.tnt]:
            mutated_mol = mutation_func(mol)
            # make sure were not getting None from mutation
            self.assertNotEqual(mutated_mol, None)

    def test_mutate_molecule(self):
        # mutate_molecule should produce 100% valid molecules (either mutated or original)
        self._run_mutation_test(self.molecular_space.mutate_molecule)

    def test_smart_atom_substitution(self):
        self._run_mutation_test(
            smart_atom_substitution,
        )

    def test_isostere_replacement(self):
        self._run_mutation_test(isostere_replacement)

    def test_functional_group_addition(self):
        self._run_mutation_test(functional_group_addition)


if __name__ == "__main__":
    unittest.main()
