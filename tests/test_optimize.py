import unittest

from src.helpers.constants import MOLECULES_SMILES
from src.Molecule import Molecule
from src.optimize import MolecularDifferentialEvolution


class TestMolecularDifferentialEvolution(unittest.TestCase):
    def setUp(self):
        # Set up a sample initial population of SMILES strings
        self.initial_smiles = [
            MOLECULES_SMILES.get("Ethanol"),  # Ethanol
            MOLECULES_SMILES.get("TNT"),  # TNT
        ]

        # Mock objective function for testing
        def mock_objective_function(molecule):
            return len(
                molecule.smiles
            )  # Fitness is based on the length of the SMILES string

        self.objective_function = mock_objective_function
        self.optimizer = MolecularDifferentialEvolution(
            self.objective_function, self.initial_smiles, max_iter=20
        )

    def test_mutation(self):
        # Test that mutation always returns a Molecule object, never None
        mutated_molecule = self.optimizer.mutate(1)
        # Assert that a Molecule object is returned
        self.assertIsNotNone(mutated_molecule, "Mutation should not return None.")
        self.assertIsInstance(
            mutated_molecule, Molecule, "Mutation should return a Molecule object."
        )

    def test_crossover(self):
        # Test the crossover function to ensure valid molecules are produced
        target = self.optimizer.molecular_space[0]
        mutant = self.optimizer.mutate(0)
        trial = self.optimizer.crossover(target, mutant)
        self.assertIsNotNone(trial, "Mutation should not return None.")
        self.assertIsInstance(
            trial, Molecule, "Crossover should return a Molecule object."
        )

    def test_selection(self):
        idx = 0
        trial = self.optimizer.molecular_space[1]
        self.optimizer.select(idx, trial)
        # Verify that the target molecule was replaced by the trial
        self.assertEqual(
            self.optimizer.molecular_space[idx].smiles,
            MOLECULES_SMILES.get("TNT"),
            "The trial molecule should replace the target molecule in the population.",
        )

        # Verify that the best fitness and best solution are updated
        self.assertEqual(
            self.optimizer.best_fitness,
            len(MOLECULES_SMILES.get("TNT")),
            "The best fitness should be updated to the trial molecule's fitness.",
        )
        self.assertEqual(
            self.optimizer.best_solution.smiles,
            MOLECULES_SMILES.get("TNT"),
            "The best solution should be updated to the trial molecule.",
        )

    def test_run_optimization(self):
        # Test if the optimization improves fitness and converges
        best_solution, best_fitness = self.optimizer.run()

        # Best solution should not be None and fitness should improve
        self.assertIsNotNone(
            best_solution, "Best solution should not be None after optimization."
        )
        self.assertGreaterEqual(
            best_fitness,
            max(len(smiles) for smiles in self.initial_smiles),
            "Best fitness should be greater or equal to intial best fitness after optimization.",
        )
        print(f"Best solution: {best_solution.smiles}, Best fitness: {best_fitness}")


if __name__ == "__main__":
    unittest.main()
