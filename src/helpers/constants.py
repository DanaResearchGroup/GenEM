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
}

ISOSTERES_LIST = [
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