from rdkit import RDLogger

from src.helpers.constants import MOLECULES_SMILES
from src.helpers.optimize_helpers import (
    normalize_mol_weight,
    normalize_ob_percentage,
    normalize_SA_score,
)
from src.optimize import MolecularDifferentialEvolution

# Set log level to suppress warnings
RDLogger.logger().setLevel(RDLogger.ERROR)


def advanced_objective_function(
    molecule, mw_weight=0.7, sascore_weight=1.2, ob_weight=0.9
):
    """
    Function to calculate the objective function with weighted properties.
    Normalizes different properties with the best score being 1.

    Arguments:
        molecule (AdvancedMolecularDifferentialEvolution): The molecule being evaluated.
        mw_weight (float): Weight for molecular weight score (default is 1.0).
        sascore_weight (float): Weight for SAScore score (default is 1.0).
        ob_weight (float): Weight for OB% score (default is 1.0).

    Returns:
        float: The molecule's fitness score based on the weighted properties.
    """

    if not molecule.properties:
        molecule.calculate_properties()

    # Normalized scores for each property
    mw_score = normalize_mol_weight(molecule.properties.get("molecular_weight", 0))
    normalized_sascore = normalize_SA_score(molecule.properties.get("sascore", 0))
    ob_score = normalize_ob_percentage(molecule.properties.get("ob_percentage", 0))

    # Apply the weights to each property score
    weighted_mw_score = mw_score * mw_weight
    weighted_sascore = normalized_sascore * sascore_weight
    weighted_ob_score = ob_score * ob_weight

    # Calculate the total weight
    total_weight = mw_weight + sascore_weight + ob_weight

    # Compute the weighted average fitness score
    fitness_score = (
        weighted_mw_score + weighted_sascore + weighted_ob_score
    ) / total_weight

    return fitness_score


# Example usage
INITIAL_SMILES = [
    MOLECULES_SMILES["Aspirin"],
    MOLECULES_SMILES["Testosterone"],
    MOLECULES_SMILES["Nicotine"],
    MOLECULES_SMILES["TNT"],
    MOLECULES_SMILES["RDX"],
    MOLECULES_SMILES["HMX"],
    MOLECULES_SMILES["PETN"],
    MOLECULES_SMILES["Tetryl"],
    # TODO: find out names (if any) and add to MOLECULES_SMILES
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "CC(=O)NC1=CC=C(O)C=C1",
    "CN(C)CCN1C2=CC=CC=C2SC3=CC=CC=C13",
    "CN(C)C(=O)C1=CC2=C(N1)C=CC3=C2C=CC4=C3C=CN4C",
]

de = MolecularDifferentialEvolution(
    advanced_objective_function, INITIAL_SMILES, crossover_prob=0.5, max_iter=1000
)
best_molecule, best_fitness = de.run()
print(f"Best molecule: {best_molecule.smiles}")
print(f"Best fitness: {best_fitness}")
print(best_molecule.properties)
