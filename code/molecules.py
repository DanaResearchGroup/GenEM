from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
import properties as p
import logging
import random
import molecule_helpers

class Molecule:
    """
    Class Molecule
    Args:
        smiles (str): SMILES string
        initial_population(list): List of initial molecules
    """
    def __init__(self, smiles, initial_population=None):
        self.smiles = smiles
        self.population = [Molecule(smiles) for smiles in initial_population] if initial_population else []
        try:
            self.mol = Chem.MolFromSmiles(smiles)
            if self.mol is None:
                raise ValueError(f"Invalid SMILES string: {smiles}")
        except Exception as e:
            print(f"Error: {e}")
            self.mol = None  # Set to None if SMILES parsing fails
        self.fingerprint = None
        self.properties = {}
        self.best_solution = None
        self.best_fitness = -float('inf')

    def calculate_fingerprint(self):
        """
        Function to calculate fingerprint
        Argument:
            self (Molecule): Molecule object
        Returns:
            Molecule's fingerprint
        """
        if self.mol is not None:
            self.fingerprint = AllChem.GetMorganFingerprintAsBitVect(self.mol, radius=2, nBits=1024)

    def calculate_properties(self):
        """
        Function to calculate properties
        Argument:
            self (Molecule): Molecule object
        Returns:
            properties (dict): Dictionary of properties
        """
        if self.mol is None:
            raise ValueError("Molecule is not valid.")
        self.properties['molecular_weight'] = Descriptors.ExactMolWt(self.mol)
        self.properties['logp'] = Descriptors.MolLogP(self.mol)
        self.properties['num_h_acceptors'] = Descriptors.NumHAcceptors(self.mol)
        self.properties['num_h_donors'] = Descriptors.NumHDonors(self.mol)
        self.properties['num_rotatable_bonds'] = Descriptors.NumRotatableBonds(self.mol)
        # Calculate OB% property
        self.properties['ob_percentage'] = p.calculate_ob_percentage(self.mol)
        # Calculate SAScore
        self.properties['sascore'] = p.calculate_sascore(self.mol)

    def update_individual(self, index, new_molecule):
        """
        Function to update individual from population
        Argument:
            index (int): index of population
            new_molecule (Molecule): new molecule
        Returns:
            new_molecule (Molecule): updated molecule
        """
        self.population[index] = new_molecule

    def calculate_population_fitness(self, objective_function):
        """
        Function to calculate population fitness using a provided objective function
        Argument:
            objective_function (function): objective function
        Returns:
            None
        """
        for mol in self.population:
            fitness = objective_function(mol)
            if fitness > self.best_fitness:
                self.best_fitness = fitness
                self.best_solution = mol

    def smart_atom_substitution(self, molecule):
        """
        Function to substitute atom with smart atom in order to increase chance of molecular validity.
        It defines allowed substitutions based on atom type and valence
        Argument:
            self(Molecule): Molecule
            molecule(Molecule): Molecule which is represented as fingerprint
        Return:
            mol(Molecule): Molecule which is represented as fingerprint
        """
        mol = Chem.RWMol(molecule.mol)
        if mol.GetNumAtoms() > 0:
            atom_idx = random.randint(0, mol.GetNumAtoms() - 1)
            atom = mol.GetAtomWithIdx(atom_idx)

            SUBSTITUTIONS = {
                'C': ['N', 'O', 'S'],
                'N': ['C', 'O', 'S'],
                'O': ['N', 'S'],
                'S': ['O', 'N']
            }

            original_atom_symbol = atom.GetSymbol()
            if original_atom_symbol in SUBSTITUTIONS:
                new_atom_symbol = random.choice(SUBSTITUTIONS[original_atom_symbol])
                new_atom = Chem.Atom(new_atom_symbol)
                new_atom.SetFormalCharge(0)  # Reset formal charge
                mol.ReplaceAtom(atom_idx, new_atom)

        return self.sanitize_and_optimize_molecule(mol)

    def isostere_replacement(self, molecule):
        """
        Function to replace similar functional groups with other functional groups.
        Argument:
        self(Molecule): Molecule
        molecule(Molecule): Molecule which is represented as fingerprint
        Return:
        mol(Molecule): Molecule which is represented as fingerprint
        """
        original_smiles = Chem.MolToSmiles(molecule.mol)
        mol = Chem.RWMol(molecule.mol)
        ISOSTERES = {
            "C(=O)O": "C(=O)NH2",  # Carboxylic acid to amide
            "C(F)(F)F": "C#N",  # Trifluoromethyl to nitrile
            "c1ccccc1": "c1ccncc1",  # Benzene to pyridine
            "C1=CC2=C(C=C1)C(=O)N2C(=O)O": "C1=CC2=C(C=C1)C(=O)N2C(=O)O",
            # Nitro group replacement for energetic materials
            "C[N+](=O)[O-]": "C1=CC=C(C(=O)O)C=C1",  # Nitro to carboxyl group
            # Nitro group replacement to fused ring system
            "[O-][N+](=O)C": "[O-][N+](=O)C1=CC=CC=C1",  # Nitro group to phenyl ring
            "[N+](=O)[O-]": "[N+](=O)[O-]C1=CC=C(C(=O)O)C=C1"  # Nitro group to substituted phenyl ring
        }

        try:
            for original, replacement in ISOSTERES.items():
                patt = Chem.MolFromSmarts(original)
                repl = Chem.MolFromSmiles(replacement)

                if patt and repl and mol.HasSubstructMatch(patt):
                    mol = AllChem.ReplaceSubstructs(mol, patt, repl, replaceAll=True)[0]
                    break
        except Exception:
            return None  # Return None if the replacement fails

        smiles = Chem.MolToSmiles(mol)

        # Check if the molecule was changed
        if smiles == original_smiles:
            return None

        # Combine fragments if necessary
        if '.' in smiles:
            try:
                combined_mol = self.combine_fragments(smiles)
                if combined_mol:
                    mol = combined_mol
            except Exception:
                return None

        return self.sanitize_and_optimize_molecule(mol)

    def functional_group_addition(self, molecule):
        """
        Function to add functional groups to molecule.
        Argument:
            self(Molecule): Molecule
            molecule(Molecule): Molecule which is represented as fingerprint
        Return:
            mol(Molecule): Molecule which is represented as fingerprint
        """
        mol = Chem.RWMol(molecule.mol)
        FUNCTIONAL_GROUPS = ["[OH]", "[NH2]", "[C](=O)[OH]", "[CH3]"]
        if mol.GetNumAtoms() > 0:
            atom_idx = random.randint(0, mol.GetNumAtoms() - 1)
            fg = Chem.MolFromSmarts(random.choice(FUNCTIONAL_GROUPS))
            if fg:
                edmol = Chem.EditableMol(mol)
                edmol.ReplaceAtom(atom_idx, fg.GetAtomWithIdx(0))
                mol = edmol.GetMol()

        return self.sanitize_and_optimize_molecule(mol)


    def sanitize_and_optimize_molecule(self, mol):
        """
        Function to sanitize the molecule after mutation to improve validity. In addition, logging has been designed
        to track invalid molecules and what mutations are causing them.
        Argument:
            self(Molecule): Molecule which is represented as fingerprint
            mol(Molecule): Molecule which is represented as fingerprint after mutation
        Return:
            Molecule(new_smiles): Molecule which is represented as fingerprint after undergoing sanitization.
        """
        try:
            Chem.SanitizeMol(mol)
            mol = Chem.AddHs(mol)
            new_smiles = Chem.MolToSmiles(mol)
            return Molecule(new_smiles)
        except Exception:
            return None