from optimize import MolecularDifferentialEvolution


def advanced_objective_function(molecule):
    """
    Function to calculate the objective function. Normalizes different properties with the best score being 1.
    Arguments:
        molecule(AdvancedMolecularDifferentialEvolution)
    Returns:
        molecule's fitness score
    """
    if not molecule.properties:
        molecule.calculate_properties()

    mw_score = 1 - abs(molecule.properties['molecular_weight'] - 300) / 300
    logp_score = 1 - abs(molecule.properties['logp'] - 2) / 5
    normalized_sascore = (10 - molecule.properties['sascore']) / 9

    return (mw_score + logp_score + normalized_sascore) / 3

# Example usage
initial_smiles = [
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

de = MolecularDifferentialEvolution(advanced_objective_function, initial_smiles, crossover_prob=0.5, max_iter=500)
best_molecule, best_fitness = de.run()
print(f"Best molecule: {best_molecule.smiles}")
print(f"Best fitness: {best_fitness}")
