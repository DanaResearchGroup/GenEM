import unittest

from rdkit import Chem

from src.helpers.constants import MOLECULES_SMILES
from src.molecule import Molecule


class TestMolecule(unittest.TestCase):

    def setUp(self):
        self.ethanol_smiles = MOLECULES_SMILES.get("Ethanol")  # Ethanol
        self.tnt_smiles = MOLECULES_SMILES.get("TNT")  # TNT
        self.invalid_smiles = "COCOCOCOCOC11##111OOO"  # invalid smile
        self.ethanol = Molecule(self.ethanol_smiles)
        self.tnt = Molecule(self.tnt_smiles)

    def test_smiles_to_mol(self):
        expected_eth_smiles = Chem.CanonSmiles(
            self.ethanol_smiles
        )  # Canonicalize SMILES
        generated_eth_smiles = Chem.CanonSmiles(Chem.MolToSmiles(self.ethanol.mol))
        expected_tnt_smiles = Chem.CanonSmiles(self.tnt_smiles)  # Canonicalize SMILES
        generated_tnt_smiles = Chem.CanonSmiles(Chem.MolToSmiles(self.tnt.mol))
        self.assertEqual(generated_eth_smiles, expected_eth_smiles)
        self.assertEqual(generated_tnt_smiles, expected_tnt_smiles)

    def test_smart_atom_substitution(self):
        for mol in [self.ethanol, self.tnt]:
            mutated_mol = Molecule.smart_atom_substitution(mol)
            # make sure were not getting None from mutation
            if mutated_mol is not None:
                print(Chem.MolToSmiles(mutated_mol.mol))
                self.assertNotEqual(mol.mol.ToBinary(), mutated_mol.mol.ToBinary())

    def test_isostere_replacement(self):
        for mol in [self.ethanol, self.tnt]:
            mutated_mol = Molecule.isostere_replacement(mol)
            # make sure were not getting None from mutation
            if mutated_mol is not None:
                print(Chem.MolToSmiles(mutated_mol.mol))
                self.assertNotEqual(mol.mol.ToBinary(), mutated_mol.mol.ToBinary())

    def test_functional_group_addition(self):
        for mol in [self.ethanol, self.tnt]:
            mutated_mol = Molecule.functional_group_addition(mol)
            # make sure were not getting None from mutation
            if mutated_mol is not None:
                print(Chem.MolToSmiles(mutated_mol.mol))
                self.assertNotEqual(mol.mol.ToBinary(), mutated_mol.mol.ToBinary())


if __name__ == "__main__":
    unittest.main()
