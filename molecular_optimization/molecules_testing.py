#!/usr/bin/env python3
# encoding: utf-8

import unittest
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from molecules import Molecule


class TestMolecule(unittest.TestCase):

    def setUp(self):
        self.ethanol_smiles = "CCO"  # Ethanol
        self.tnt_smiles = "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]"  # TNT
        self.invalid_smiles = "COCOCOCOCOC11##111OOO" #invalid smile
        self.ethanol = Molecule(self.ethanol_smiles)
        self.tnt = Molecule(self.tnt_smiles)
        self.invalid_molecule = Molecule(self.invalid_smiles)

    def test_smiles_to_mol(self):
        expected_eth_smiles = Chem.CanonSmiles(self.ethanol_smiles)  # Canonicalize SMILES
        generated_eth_smiles = Chem.CanonSmiles(Chem.MolToSmiles(self.ethanol.mol))
        expected_tnt_smiles = Chem.CanonSmiles(self.tnt_smiles)  # Canonicalize SMILES
        generated_tnt_smiles = Chem.CanonSmiles(Chem.MolToSmiles(self.tnt.mol))
        self.assertEqual(generated_eth_smiles, expected_eth_smiles)
        self.assertEqual(generated_tnt_smiles, expected_tnt_smiles)

    def test_fingerprint(self):
        self.ethanol.calculate_fingerprint()
        self.assertIsNotNone(self.ethanol.fingerprint)
        self.assertEqual(len(self.ethanol.fingerprint), 1024)  # Default nBits=1024

    def test_calculate_properties(self):
        # Test property calculation
        self.ethanol.calculate_properties()
        props_eth = self.ethanol.properties
        self.assertIn('molecular_weight', props_eth)
        self.assertIn('logp', props_eth)
        self.assertIn('num_h_acceptors', props_eth)
        self.assertIn('num_h_donors', props_eth)
        self.assertIn('num_rotatable_bonds', props_eth)
        self.assertAlmostEqual(props_eth['molecular_weight'], 46.04186, places=4)  # Ethanol MW
        self.assertAlmostEqual(props_eth['logp'], -0.0014, places=4)  # LogP of ethanol
        self.assertGreaterEqual(props_eth['sascore'], 1)
        self.assertLessEqual(props_eth['sascore'],10)

        self.tnt.calculate_properties()
        # Check if properties are calculated correctly for TNT
        self.assertAlmostEqual(self.tnt.properties['molecular_weight'], 227.0, places=1)
        self.assertAlmostEqual(self.tnt.properties['logp'], 1.7, places=1)  # Approximate value
        self.assertEqual(self.tnt.properties['num_h_acceptors'], 6)
        self.assertEqual(self.tnt.properties['num_h_donors'], 0)
        self.assertEqual(self.tnt.properties['num_rotatable_bonds'], 3)

    def test_generate_3d_conformer(self):
        self.ethanol.generate_3d_conformer()
        self.tnt.generate_3d_conformer()
        self.assertIsNotNone(self.ethanol.conformer)
        self.assertIsNotNone(self.tnt.conformer)

        # Check that the invalid molecule cannot generate a conformer
        self.assertIsNone(self.invalid_molecule.mol)

    def test_calculate_energy(self):
        eth_energy = self.ethanol.calculate_energy()
        tnt_energy = self.tnt.calculate_energy()
        self.assertIsInstance(eth_energy, float)  # Ensure energy is a float
        self.assertGreater(eth_energy, -1000)  # Ensure it's not a very large negative value
        self.assertIsInstance(tnt_energy, float)  # Ensure energy is a float
        self.assertGreater(tnt_energy, -1000)

    def test_invalid_molecule_energy(self):
        # Test energy calculation for invalid molecule (None case)
        invalid_molecule = Molecule("invalid_smiles")
        energy = invalid_molecule.calculate_energy()
        self.assertEqual(energy, float('inf'))  # Energy should be inf if molecule invalid


if __name__ == '__main__':
    unittest.main()
