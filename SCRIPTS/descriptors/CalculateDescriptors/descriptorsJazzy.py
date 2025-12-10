import pandas
from rdkit import Chem
from rdkit.Chem import Crippen
from jazzy.api import deltag_from_smiles
from jazzy.api import molecular_vector_from_smiles
import warnings


def calculateDescriptorsJazzy(smiles):
    warnings.filterwarnings("ignore")
    columnNames = {"sdc":"CHds","sdx":"XHds","sa":"HBAs"}
    energyColumnNames = {"dga":"HydA","dgp":"HydP","dgtot":"Hyd"}
    deltag_sol_oct = []
    v_dict = {}
    descriptors = {}
    for smile in smiles:
        try:
            canonSmile = Chem.CanonSmiles(smile)
            molecule = Chem.MolFromSmiles(canonSmile)
            v = molecular_vector_from_smiles(canonSmile)
            for key in columnNames:
                if columnNames[key] in v_dict:
                    v_dict[columnNames[key]].append(v[key])
                else:
                    v_dict[columnNames[key]] = [v[key]]
            for key in energyColumnNames:
                if energyColumnNames[key] in v_dict:
                    v_dict[energyColumnNames[key]+"0"].append(v[key])
                    v_dict[energyColumnNames[key]].append(v[key]/4.184)
                else:
                    v_dict[energyColumnNames[key]+"0"] = [v[key]]
                    v_dict[energyColumnNames[key]] = [v[key]/4.184]
            logP = Chem.Crippen.MolLogP(molecule)
            d = deltag_from_smiles(canonSmile)/4.184
            dg_sol_oct = -logP*1.36+d
            deltag_sol_oct.append(dg_sol_oct)
        except:
            print(f"Jazzy: Failed to calculate descriptors for {smile}")
            deltag_sol_oct.append(" ")
            for key in columnNames:
                if columnNames[key] in v_dict:
                    v_dict[columnNames[key]].append(" ")
                else:
                    v_dict[columnNames[key]] = [" "]
            for key in energyColumnNames:
                if energyColumnNames[key] in v_dict:
                    v_dict[energyColumnNames[key]+"0"].append(" ")
                    v_dict[energyColumnNames[key]].append(" ")
                else:
                    v_dict[energyColumnNames[key]+"0"] = [" "]
                    v_dict[energyColumnNames[key]] = [" "]
    for key in v_dict:
        descriptors[key] = v_dict[key]
    descriptors["delta_g_sol_oct"] = deltag_sol_oct
    return pandas.DataFrame(descriptors)

