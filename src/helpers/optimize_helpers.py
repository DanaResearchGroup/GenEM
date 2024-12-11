import random

from src.molecule import Molecule


def mutate_molecule(molecule: Molecule) -> Molecule:
    """
    Function to mutate a molecule. It works by going through the list of mutations, randomly picking one. If the
    mutation results in a valid molecule, it returns the mutated molecule. Otherwise, it will randomly pick another
    mutation, and repeat the process. If the mutation doesn't produce a valid molecule after 3 times,
    it returns the original molecule.
    Argument:
        molecule(Molecule): Molecule which is represented as fingerprint
    Return:
        mutated_molecule(Molecule): Mutated molecule
    """
    mutation_strategies = [
        molecule.smart_atom_substitution,
        molecule.backbone_replacement,
        molecule.functional_group_addition,
    ]
    for _ in range(3):
        strategy = random.choice(mutation_strategies)
        mutated_molecule = strategy(molecule)
        if mutated_molecule and mutated_molecule.rdkit_mol:
            return mutated_molecule
    return molecule


def normalize_property(score, ideal, max_range):
    return max(0, 1 - abs(score - ideal) / max_range)
