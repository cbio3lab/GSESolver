from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
import os
import pandas as pd
from openbabel import pybel

def smiles2sdf(smiles, output_file):
    if smiles2sdf_rdkit(smiles, output_file):
        print(f"Failed to create sdf file with RDKit for {smiles}, trying with Obabel")
        return True
    elif smiles2sdf_obabel(smiles, output_file):
        return True
    else: 
        return False


def smiles2sdf_rdkit(smiles, output_file):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return
    mol = Chem.AddHs(mol)
    count = 0
    max_count = 20
    is_embedding_succesful = False
    while count < max_count:
        result = AllChem.EmbedMolecule(mol, randomSeed=count, maxAttempts=5000)
        count += 1
        if result == 0:
            is_embedding_succesful = True
            break
    if is_embedding_succesful:
        AllChem.UFFOptimizeMolecule(mol)
        writer = Chem.SDWriter(output_file)
        writer.write(mol)
        writer.close()
        return True
    else:
        return False

def smiles2sdf_obabel(smiles, output_file):
    try:
        mol = pybel.readstring("smi",smiles)
        mol.addh()
        mol.make3D()
        mol.write("sdf",output_file,overwrite=True)
        return True
    except:
        return False



def calculateConjugatedDescriptors(smiles_list):
    for i, smiles in enumerate(smiles_list):
        output_file = f'molecule_{i+1}.sdf'
        is_sdf_creation_successful = smiles2sdf(smiles, output_file)
        if not is_sdf_creation_successful:
            print(f"Failed at creating sdf for {smiles}")
    subprocess.call("Rscript conjugaR.R", shell=True)
    files = os.listdir()
    sdf_files = [file for file in files if file.split(".")[-1] == "sdf"]
    descriptors_file_path = 'conjugated_descriptors.csv'
    descriptors = pd.read_csv(descriptors_file_path)
    descriptors.drop(columns=[descriptors.columns[0]],axis=1,inplace=True)
    for file in sdf_files:
        os.remove(file)
    for i in range(len(smiles_list)):
        try:
            file = f"molecule_{i+1}_conjugate_counts.csv"
            os.remove(file)
        except:
            pass
    os.remove(descriptors_file_path)
    return descriptors

