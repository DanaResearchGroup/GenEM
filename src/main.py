from rdkit import RDLogger
from rdkit.Chem import Descriptors

from src.helpers.constants import MOLECULES_SMILES
from src.helpers.optimize_helpers import normalize_property
from src.optimize import MolecularDifferentialEvolution

# Set log level to suppress warnings
RDLogger.logger().setLevel(RDLogger.ERROR)


def advanced_objective_function(mol):
    """
    Function to calculate the objective function with weighted properties.
    Normalizes different properties with the best score being 1.

    Arguments:
        mol (Molecule): The molecule being evaluated.
    Returns:
        float: The molecule's fitness score based on the weighted properties.
    """

    # Calculate 'fake' properties to imitate optimization
    mw_score = Descriptors.ExactMolWt(mol.rdkit_mol)
    logp_score = Descriptors.MolLogP(mol.rdkit_mol)
    rot_bonds_score = Descriptors.NumRotatableBonds(mol.rdkit_mol)
    h_acc_score = Descriptors.NumHAcceptors(mol.rdkit_mol)
    rads_score = Descriptors.NumRadicalElectrons(mol.rdkit_mol)
    # Normalize the properties between 1 and 0 where 1 is best
    # Example normalization
    normalized_mw_score = normalize_property(mw_score, ideal=250, max_range=250)
    normalized_logp_score = normalize_property(
        logp_score, ideal=0, max_range=5
    )  # Assuming max LogP deviation is 5
    normalized_rot_bonds_score = normalize_property(
        rot_bonds_score, ideal=0, max_range=10
    )
    normalized_h_acc_score = normalize_property(h_acc_score, ideal=1, max_range=5)
    normalized_rads_score = normalize_property(rads_score, ideal=0, max_range=5)

    # Compute the average fitness score
    fitness_score = (
        normalized_mw_score
        + normalized_rads_score
        + normalized_h_acc_score
        + normalized_rot_bonds_score
        + normalized_logp_score
    ) / 5

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
    MOLECULES_SMILES["Caffeine"],
    MOLECULES_SMILES["Ibuprofen"],
    MOLECULES_SMILES["Acetaminophen"],
    MOLECULES_SMILES["Promethazine"],
    MOLECULES_SMILES["Harmine"]
]

de = MolecularDifferentialEvolution(
    advanced_objective_function, INITIAL_SMILES, crossover_prob=0.5, max_iter=50
)
best_molecule, best_fitness = de.run()
print("\nBest Final Fitness:")
print(f"SMILES={best_molecule.smiles}, Fitness={best_fitness}")
