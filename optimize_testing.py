#!/usr/bin/env python3
# encoding: utf-8

import unittest
from molecules import Molecule
from optimize import AdvancedPopulation, MolecularDifferentialEvolution
from rdkit import Chem


class TestAdvancedPopulation(unittest.TestCase):

    def setUp(self):
        # Set up a sample initial population of SMILES strings
        self.initial_smiles = [
            "CCO",  # Ethanol
            "CC(C)O",  # Isopropanol
            "CCCC",  # Butane
            "CCN",  # Ethylamine
        ]

        # Create the AdvancedPopulation instance
        self.population = AdvancedPopulation(self.initial_smiles)

        # Mock objective function for testing
        def mock_objective_function(molecule):
            return len(molecule.smiles)  # Fitness is based on the length of the SMILES string

        self.objective_function = mock_objective_function

    def test_get_individuals(self):
        # Test if we can retrieve individuals correctly
        individual = self.population.get_individuals(0)
        self.assertEqual(individual.smiles, self.initial_smiles[0],
                         "The SMILES string of the individual should match the initial SMILES at index 0.")

        individual = self.population.get_individuals(2)
        self.assertEqual(individual.smiles, self.initial_smiles[2],
                         "The SMILES string of the individual should match the initial SMILES at index 2.")

    def test_update_individual(self):
        # Test updating an individual in the population
        new_molecule = Molecule("CO")  # Methanol

        self.population.update_individual(0, new_molecule)

        updated_individual = self.population.get_individuals(0)
        self.assertEqual(updated_individual.smiles, "CO", "The SMILES string of the updated individual should be 'CO'.")
        self.assertNotEqual(updated_individual.smiles, self.initial_smiles[0],
                            "The updated individual should no longer have the original SMILES.")

    def test_calculate_population_fitness(self):
        # Test if the fitness is calculated correctly
        self.population.calculate_population_fitness(self.objective_function)

        # The best fitness should correspond to the molecule with the longest SMILES string
        expected_best_fitness = max(len(smiles) for smiles in self.initial_smiles)
        expected_best_solution_smiles = max(self.initial_smiles, key=len)

        self.assertEqual(self.population.best_fitness, expected_best_fitness,
                         "Best fitness should correspond to the molecule with the longest SMILES string.")
        self.assertEqual(self.population.best_solution.smiles, expected_best_solution_smiles,
                         "Best solution should correspond to the molecule with the longest SMILES string.")


class TestMolecularDifferentialEvolution(unittest.TestCase):

    def setUp(self):
        # Set up a sample initial population of SMILES strings
        self.initial_population = [
            "CCO",  # Ethanol
            "CC(C)O",  # Isopropanol
            "CCCC",  # Butane
            "CCN",  # Ethylamine
        ]

        # Dummy objective function for testing
        def mock_objective_function(molecule):
            return len(molecule.smiles)  # Simple objective based on SMILES length

        self.objective_function = mock_objective_function
        self.optimizer = MolecularDifferentialEvolution(self.objective_function, self.initial_population, max_iter=20)

    def test_mutation(self):
        # Test the mutation function to ensure it returns a valid molecule
        mutated_molecule = self.optimizer.mutate(0)

        self.assertIsInstance(mutated_molecule, Molecule, "Mutation should return a Molecule object.")

    def test_crossover(self):
        # Test the crossover function to ensure valid molecules are produced
        target = self.optimizer.advanced_molecular_space.population[0]
        mutant = self.optimizer.mutate(0)

        trial_molecule = self.optimizer.crossover(target, mutant)

        self.assertIsInstance(trial_molecule, Molecule, "Crossover should return a Molecule object.")


    def test_selection(self):
        # Test the selection function to ensure the fitter molecule is selected
        target = self.optimizer.advanced_molecular_space.population[0]
        mutant = self.optimizer.mutate(0)
        trial = self.optimizer.crossover(target, mutant)

        # Manually evaluate fitness of target and trial
        target_fitness = self.optimizer.objective_function(target)
        trial_fitness = self.optimizer.objective_function(trial)

        # Perform selection
        self.optimizer.select(0, trial)

        if trial_fitness > target_fitness:
            self.assertEqual(self.optimizer.advanced_molecular_space.population[0].smiles, trial.smiles,
                             "The trial molecule should replace the target molecule if it has better fitness.")
        else:
            self.assertEqual(self.optimizer.advanced_molecular_space.population[0].smiles, target.smiles,
                             "The target molecule should remain if it has better fitness.")

    def test_run_optimization(self):
        # Test if the optimization improves fitness and converges
        best_solution, best_fitness = self.optimizer.run()

        # Best solution should not be None and fitness should improve
        self.assertIsNotNone(best_solution, "Best solution should not be None after optimization.")
        self.assertGreaterEqual(best_fitness, max(len(smiles) for smiles in self.initial_population), "Best fitness should be greater or equal to intial best fitness after optimization.")
        print(f"Best solution: {best_solution.smiles}, Best fitness: {best_fitness}")


if __name__ == '__main__':
    unittest.main()
