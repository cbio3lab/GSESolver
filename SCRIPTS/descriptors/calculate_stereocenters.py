import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

def calculate_stereocenters(smiles_input_list: list) -> list:
    """
    Calculates the 'NumAtomStereoCenters_rdkit' descriptor for a list of SMILES strings.

    Args:
        smiles_input_list (list): A list of SMILES strings.

    Returns:
        list: A list named 'NumAtomStereoCenters_rdkit' containing the descriptor values.
              Returns None for invalid SMILES entries.
    """
    NumAtomStereoCenters_rdkit = []
    for smiles in smiles_input_list:
        mol = Chem.MolFromSmiles(smiles) if pd.notnull(smiles) else None
        if mol is not None:
            num_stereo_centers = Chem.Fragments.NumAtomStereoCenters(mol)
            NumAtomStereoCenters_rdkit.append(num_stereo_centers)
        else:
            NumAtomStereoCenters_rdkit.append(None) # Append None for invalid SMILES
    return NumAtomStereoCenters_rdkit
