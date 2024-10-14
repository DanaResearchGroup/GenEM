#!/usr/bin/env python3
# encoding: utf-8

import unittest
from mol_space import MolecularSpace
from rdkit import Chem

class TestMolecularSpace(unittest.TestCase):
    def setUp(self):
        self.molecular_space = MolecularSpace([
            'CC1=C(C(=CC(=C1)[N+](=O)[O-]))[N+](=O)[O-]',   # TNT
            'C1N2C(=O)N(C(=O)N1C2=O)C',   # RDX
            'C1C(N2C(=O)N(C(=O)N1C2=O)C)',    # HMX
            'CC(CO[N+](=O)[O-])([N+](=O)[O-])CO[N+](=O)[O-]',    # PETN
            'C1=CC(=C(C(=C1)[N+](=O)[O-])C([N+](=O)[O-])=O)'])

    def test_calculate_distance(self):
        mol1 = self.molecular_space.population[0]
        mol2 = self.molecular_space.population[1]
        distance = self.molecular_space.calculate_distance(mol1, mol2)
        self.assertGreaterEqual(distance, 0)

    def test_get_nearest_neighbors(self):
        mol = self.molecular_space.population[0]
        neighbors = self.molecular_space.get_nearest_neighbors(mol)
        self.assertEqual(len(neighbors), 4)

    def _run_mutation_test(self, mutation_func, func_name, expect_none=True):
        for i, mol in enumerate(self.molecular_space.population):
            valid_count = 0
            none_count = 0
            print(f"Testing {func_name} on molecule {i + 1}: {Chem.MolToSmiles(mol.mol)}")

            for _ in range(10):
                mutated_mol = mutation_func(mol)

                # Ensure mutated molecule is not None and is different from original
                if mutated_mol is not None and Chem.MolToSmiles(mutated_mol.mol) != Chem.MolToSmiles(mol.mol):
                    valid_count += 1
                elif mutated_mol is None:
                    none_count += 1

                print(f"Mutated mol: {Chem.MolToSmiles(mutated_mol.mol) if mutated_mol else 'None'}")

            # For individual mutation functions, we accept some None values. It is important that the mutations
            # produce less original molecules as output than both None and valid mutations. This allows for a higher
            # chance to eventually produce a valid mutant that expands exploration.
            if expect_none:
                self.assertGreaterEqual(valid_count + none_count, 5,
                                        f"{func_name} failed on molecule {i + 1}: "
                                        f"returned original molecule more than 50% of the time")
            # For mutate_molecule, we expect 0 None results
            else:
                self.assertEqual(none_count, 0, f"{func_name} returned None on molecule {i + 1}: mutation failed.")
                self.assertGreater(valid_count, 0,
                                   f"{func_name} failed on molecule {i + 1}: mutation did not change the molecule.")

    def test_mutate_molecule(self):
        # mutate_molecule should produce 100% valid molecules (either mutated or original)
        self._run_mutation_test(self.molecular_space.mutate_molecule, "mutate_molecule", expect_none=False)

    def test_smart_atom_substitution(self):
        self._run_mutation_test(self.molecular_space.smart_atom_substitution, "smart_atom_substitution")

    def test_isostere_replacement(self):
        self._run_mutation_test(self.molecular_space.isostere_replacement, "isostere_replacement")

    def test_functional_group_addition(self):
        self._run_mutation_test(self.molecular_space.functional_group_addition, "functional_group_addition")


if __name__ == '__main__':
    unittest.main()
