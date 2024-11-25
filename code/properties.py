import re

import sascorer as sa
from rdkit.Chem import Descriptors, rdMolDescriptors


def calculate_ob_percentage(mol):
    """
    Function to calculate Oxygen Bonds percentage (OB%)
    Argument:
        self (Molecule): Molecule object
    Returns:
        ob_percentage (float): Oxygen bond percentage
    """
    formula = rdMolDescriptors.CalcMolFormula(mol)

    # Regular expression to match elements and their counts
    element_counts = re.findall(r"([A-Z][a-z]*)(\d*)", formula)

    # Initialize counts for C, H, O, and M
    x = 0  # Carbon count
    y = 0  # Hydrogen count
    z = 0  # Oxygen count
    m = 0  # Metal oxide count (assuming 0 at the moment)

    # Loop through the element counts and assign values to C, H, O
    for element, count in element_counts:
        count = int(count) if count else 1  # Default to 1 if no subscript is present
        if element == "C":
            x = count
        elif element == "H":
            y = count
        elif element == "O":
            z = count

    mw = Descriptors.ExactMolWt(mol)
    ob_percentage = (-1600 / mw) * (2 * x + y / 2 + m - z)
    return ob_percentage


def calculate_sascore(mol):
    try:
        sascore = sa.calculateScore(mol)
    except Exception as e:
        print(f"Error calculating SAScore: {e}")
        sascore = None
    return sascore
