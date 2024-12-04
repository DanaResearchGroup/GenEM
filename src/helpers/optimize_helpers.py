import random

from Molecule import Molecule


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
        molecule.isostere_replacement,
        molecule.functional_group_addition,
    ]
    for _ in range(3):
        strategy = random.choice(mutation_strategies)
        mutated_molecule = strategy(molecule)
        if mutated_molecule and mutated_molecule.mol:
            return mutated_molecule
    return molecule


def normalize_mol_weight(mol_weight: int):
    """
    Function to normalize molecular weight.

    Ideal molecular weight is 250.
    """
    return 1 - abs(mol_weight - 250) / 250


def normalize_SA_score(sa_score: int):
    """
    Function to normalize molecular weight.

    Ideal SAScore is 1, worst is 10.
    """
    return (10 - sa_score) / 9


def normalize_ob_percentage(ob_percentage: int):
    """
    Function to normalize molecular weight.

    Ideal OB% is 0%.

    # Typical OB% for EMs:
    # TNT(Trinitrotoluene): around - 74 % (fuel - rich)
    # RDX(Cyclotrimethylenetrinitramine): around - 22 %
    # HMX(Cyclotetramethylenetetranitramine): around - 21 %
    # PETN(Pentaerythritol tetranitrate): around 0 % (near - optimal balance)
    """
    return 1 - abs(ob_percentage) / 100
