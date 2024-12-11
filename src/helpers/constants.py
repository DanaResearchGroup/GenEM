MOLECULES_SMILES = {
    "Ethanol": "CCO",
    "TNT": "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]",
    "RDX": "C1N2C(=O)N(C(=O)N1C2=O)C",
    "HMX": "C1C(N2C(=O)N(C(=O)N1C2=O)C)",
    "PETN": "CC(CO[N+](=O)[O-])([N+](=O)[O-])CO[N+](=O)[O-]",
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Testosterone": "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C",
    "Nicotine": "CN1CCC[C@H]1C2=CN=CC=C2",
    "Tetryl": "C1=CC(=C(C(=C1)[N+](=O)[O-])C([N+](=O)[O-])=O)",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "Acetaminophen": "CC(=O)NC1=CC=C(O)C=C1",
    "Promethazine": "CN(C)CCN1C2=CC=CC=C2SC3=CC=CC=C13",
    "Harmine": "CN(C)C(=O)C1=CC2=C(N1)C=CC3=C2C=CC4=C3C=CN4C",
}

BACKBONE_LIST = [
    "C(=O)O",
    "C(=O)N",
    "C(F)(F)F",
    "C#N",
    "c1ccccc1",
    "c1ccncc1",
    "C1=CC2=C(C=C1)C(=O)N2C(=O)O",
    "C1=CC2=C(C=C1)C(=O)N2C(=O)O",
    "C[N+](=O)[O-]",
    "CC1=CC=C(C(=O)O)C=C1",
    "[O-][N+](=O)C",
    "[O-][N+](=O)C1=CC=CC=C1",
    "[N+](=O)[O-]",
    # "[N+](=O)[O-]C1=CC=C(C(=O)O)C=C1",
]

FUNCTIONAL_GROUPS = ["[OH]", "[NH2]", "[C](=O)[OH]", "[CH3]", "NN=O", "C=C", "CC", "NN", "N=N", "N", "O", "C"]

SUBSTITUTIONS = ["C", "N", "O", "S", "H"]