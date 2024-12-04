import numpy as np
from rdkit import Chem

import helpers.optimize_helpers as oh
from Molecule import Molecule


class MolecularDifferentialEvolution:
    """
    Class MolecularDifferentialEvolution
    Argument:
        objective_function (function): objective function
        initial_smiles (list): initial smiles will represent the initial population for optimization
        crossover_prob (float): crossover probability
        max_iter (int): maximum number of iterations
    """

    def __init__(
        self, objective_function, initial_population, crossover_prob=0.7, max_iter=1000
    ):
        self.objective_function = objective_function
        self.population_size = len(initial_population)
        self.crossover_prob = crossover_prob
        self.max_iter = max_iter
        self.molecular_space = [Molecule(smiles) for smiles in initial_population]
        self.best_solution = None
        self.best_fitness = -float(
            "inf"
        )  # Start with negative infinity for maximization goal
        # Print initial fitness values for each molecule
        self.print_initial_fitness()

    def print_initial_fitness(self):
        """
        Function for printing initial fitness values of population
        Returns:
            prints initial fitness values of population
        """
        print("Initial fitness values of starting molecules:")
        for i, molecule in enumerate(self.molecular_space):
            fitness = self.objective_function(molecule)
            print(f"Molecule {i+1}: SMILES={molecule.smiles}, Fitness={fitness}")

    def mutate(self, idx):
        """
        Function to mutate molecule
        Argument:
            idx (int): index of molecule
        Returns:
            mutated molecule if mutation was successful, otherwise target (original molecule)
        """
        target = self.molecular_space[idx]
        mutant = oh.mutate_molecule(target)

        if mutant is None:
            return target  # Fallback to target molecule if mutation fails

        return mutant

    def crossover(self, target, mutant):
        """
        Function for molecule crossover
        Argument:
            target (AdvancedMolecule): original molecule
            mutant (AdvancedMolecule): mutated molecule
        Returns:
            trial molecule if crossover was successful, otherwise target (original molecule)
        """
        try:
            target_smiles = target.smiles
            mutant_smiles = mutant.smiles

            target_mol = Chem.MolFromSmiles(target_smiles)
            mutant_mol = Chem.MolFromSmiles(mutant_smiles)

            if target_mol is None or mutant_mol is None:
                return target  # Fallback to target if conversion fails

            crossover_mask = (
                np.random.rand(target_mol.GetNumAtoms()) < self.crossover_prob
            )
            trial_mol = Chem.RWMol(target_mol)

            for i in range(trial_mol.GetNumAtoms()):
                if crossover_mask[i]:
                    trial_mol.ReplaceAtom(i, mutant_mol.GetAtomWithIdx(i))

            Chem.SanitizeMol(trial_mol)  # Try to sanitize
            trial_smiles = Chem.MolToSmiles(trial_mol)
            return Molecule(trial_smiles)
        except Exception:
            return target  # Return target if trial molecule is invalid

    def select(self, idx, trial):
        """
        Function to select a molecule by fitness score of objective function
        Argument:
            idx (int): index of molecule
            trial (AdvancedMolecule): trial molecule, which was created during crossover by combined mutated molecule
            and original molecule
        Returns:
            trial molecule if its fitness was better than the target molecule, otherwise target (original molecule)
        """
        target = self.molecular_space[idx]
        target_fitness = self.objective_function(target)
        trial_fitness = self.objective_function(trial)

        if trial_fitness > target_fitness:
            self.molecular_space[idx] = trial

        if trial_fitness > self.best_fitness:
            self.best_fitness = trial_fitness
            self.best_solution = trial

    def run(self):
        """
        Function to run optimization
        Returns:
            best_solution (AdvancedMolecule)
            best_fitness (float)
            count (int): number of unsuccessful crossovers
        """
        for _ in range(self.max_iter):
            for i in range(self.population_size):
                target = self.molecular_space[i]
                mutant = self.mutate(i)
                self.crossover(target, mutant)
                self.select(i, mutant)

        return self.best_solution, self.best_fitness
