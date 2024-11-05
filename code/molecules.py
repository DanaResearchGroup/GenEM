from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
import re
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

    def calculate_ob_percentage(self):
        """
        Function to calculate Oxygen Bonds percentage (OB%)
        Argument:
            self (Molecule): Molecule object
        Returns:
            ob_percentage (float): Oxygen bond percentage
        """
        formula = rdMolDescriptors.CalcMolFormula(self.mol)

        # Regular expression to match elements and their counts
        element_counts = re.findall(r'([A-Z][a-z]*)(\d*)', formula)

        # Initialize counts for C, H, O, and M
        carbon_count = 0  
        hyrdro_count = 0  
        oxy_count = 0
        metal_count = 0  #(assuming 0 at the moment)

        # Loop through the element counts and assign values to C, H, O
        for element, count in element_counts:
            count = int(count) if count else 1  # Default to 1 if no subscript is present
            if element == 'C':
                carbon_count = count
            elif element == 'H':
                hyrdro_count = count
            elif element == 'O':
                oxy_count = count

        mol_weight = self.properties['molecular_weight']  # molecular weight
        ob_percentage = (-1600 / mol_weight) * (2 * carbon_count + hyrdro_count / 2 + metal_count - oxy_count)
        return ob_percentage

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
        self.properties['ob_percentage'] = self.calculate_ob_percentage()
        # Calculate SAScore
        try:
            self.properties['sascore'] = sa.calculateScore(self.mol)
        except Exception as e:
            print(f"Error calculating SAScore: {e}")
            self.properties['sascore'] = None
