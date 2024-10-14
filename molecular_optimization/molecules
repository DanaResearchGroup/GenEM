from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import sascorer as sa

class Molecule:
    """
    Class Molecule
    Args:
        smiles (str): SMILES string
    """

    def __init__(self, smiles):
        self.smiles = smiles
        try:
            self.mol = Chem.MolFromSmiles(smiles)
            if self.mol is None:
                raise ValueError(f"Invalid SMILES string: {smiles}")
        except Exception as e:
            print(f"Error: {e}")
            self.mol = None  # Set to None if SMILES parsing fails
        self.fingerprint = None
        self.properties = {}

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
        # Calculate SAScore
        try:
            self.properties['sascore'] = sa.calculateScore(self.mol)
        except Exception as e:
            print(f"Error calculating SAScore: {e}")
            self.properties['sascore'] = None
    def generate_3d_conformer(self):
        """
        Function to generate 3D conformer
        Argument:
            self (Molecule): Molecule object
        """
        if self.mol is None:
            return None
        mol_with_h = Chem.AddHs(self.mol)  # Add explicit hydrogens
        AllChem.EmbedMolecule(mol_with_h, randomSeed=42)  # Generate 3D coordinates
        AllChem.MMFFOptimizeMolecule(mol_with_h)  # Optimize the 3D structure
        self.conformer = mol_with_h.GetConformer()
        self.mol = mol_with_h

    def calculate_energy(self):
        """
        Function to calculate energy
        Argument:
            self (Molecule): Molecule object
        Returns:
            energy (float): Energy of the molecule
        """
        if self.mol is None:
            return float('inf')  # Return a high value if molecule is invalid

        if self.mol.GetNumConformers() == 0:
            self.generate_3d_conformer()  # Generate a 3D conformer if no conformers exist

        if self.mol.GetNumConformers() == 0:
            return float('inf')  # If conformer generation failed, return high energy for invalid molecule

        # Set up MMFF force field
        ff = AllChem.MMFFGetMoleculeForceField(self.mol, AllChem.MMFFGetMoleculeProperties(self.mol))
        if not ff:
            return float('inf')  # Return a high value if the force field setup fails

        # Calculate energy
        energy = ff.CalcEnergy()
        return energy
