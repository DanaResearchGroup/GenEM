from optimize import MolecularDifferentialEvolution

def advanced_objective_function(molecule, mw_weight=0.7, sascore_weight=1.2, ob_weight=0.9):
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
    mw_score = 1 - abs(molecule.properties['molecular_weight'] - 250) / 250  # Ideal molecular weight is 250
    normalized_sascore = (10 - molecule.properties['sascore']) / 9  # Ideal SAScore is 1, worst is 10
    ob_score = 1 - abs(molecule.properties['ob_percentage']) / 100  # Ideal OB% is 0%
    #Typical OB% for EMs:
    #TNT(Trinitrotoluene): around - 74 % (fuel - rich)
    #RDX(Cyclotrimethylenetrinitramine): around - 22 %
    #HMX(Cyclotetramethylenetetranitramine): around - 21 %
    #PETN(Pentaerythritol tetranitrate): around 0 % (near - optimal balance)

    # Apply the weights to each property score
    weighted_mw_score = mw_score * mw_weight
    weighted_sascore = normalized_sascore * sascore_weight
    weighted_ob_score = ob_score * ob_weight

    # Calculate the total weight
    total_weight = mw_weight + sascore_weight + ob_weight

    # Compute the weighted average fitness score
    fitness_score = (weighted_mw_score + weighted_sascore + weighted_ob_score) / total_weight

    return fitness_score

# Example usage
INITIAL_SMILES = [
    "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
    "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C",  # Testosterone
    "CN1CCC[C@H]1C2=CN=CC=C2",  # Nicotine
    'CC1=C(C(=CC(=C1)[N+](=O)[O-]))[N+](=O)[O-]',   #TNT
    'C1N2C(=O)N(C(=O)N1C2=O)C',   #RDX
    'C1C(N2C(=O)N(C(=O)N1C2=O)C)',    #HMX
    'CC(CO[N+](=O)[O-])([N+](=O)[O-])CO[N+](=O)[O-]',    #PETN
    'C1=CC(=C(C(=C1)[N+](=O)[O-])C([N+](=O)[O-])=O)',    #Tetryl
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "CC(=O)NC1=CC=C(O)C=C1",
    "CN(C)CCN1C2=CC=CC=C2SC3=CC=CC=C13",
    "CN(C)C(=O)C1=CC2=C(N1)C=CC3=C2C=CC4=C3C=CN4C"
]

de = MolecularDifferentialEvolution(advanced_objective_function, INITIAL_SMILES, crossover_prob=0.5, max_iter=1000)
best_molecule, best_fitness = de.run()
print(f"Best molecule: {best_molecule.smiles}")
print(f"Best fitness: {best_fitness}")
print(best_molecule.properties)
