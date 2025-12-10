import pandas as pd
import rpy2.robjects.packages as rpackages
from rpy2.robjects import pandas2ri
import rpy2.robjects as ro
from descriptorsRDKit import calculateDescriptorsRDKit
from descriptorsJazzy import calculateDescriptorsJazzy
from openbabel import pybel
from rdkit.Chem import CanonSmiles
from parameterReader import ParameterReader
from descriptorsConjugated import calculateConjugatedDescriptors

pandas2ri.activate()

from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
import logging
rpy2_logger.setLevel(logging.ERROR)


def calculateDescriptorsCDK(smiles):
    rcdk = rpackages.importr("rcdk")
    base = rpackages.importr("base")
    descCategories = rcdk.get_desc_categories()
    getDescNames = rcdk.get_desc_names

    descNames = base.unique(base.unlist(base.sapply(descCategories,getDescNames)))
    descNames = [name for name in descNames if name != "org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor"]
    descNames = ro.StrVector(descNames)

    mols = rcdk.parse_smiles(ro.StrVector(["C"]))
    descriptorsCDK = rcdk.eval_desc(mols,descNames)
    columns = ro.conversion.rpy2py(descriptorsCDK).columns
    descCDK = pd.DataFrame()
    i = 0
    for smi in smiles:
        try:
            mols = rcdk.parse_smiles(ro.StrVector([smi]))
            descriptorsCDK = rcdk.eval_desc(mols,descNames)

            descriptorsCDK = ro.conversion.rpy2py(descriptorsCDK)
            descCDK = pd.concat([descCDK,descriptorsCDK],ignore_index=True)
        except: 
            descriptorsCDK = pd.DataFrame(pd.NA,columns=columns,index=[i])
            descCDK = pd.concat([descCDK,descriptorsCDK],ignore_index=True)
            print(f"CDK: Failed to calculate descriptors for {smi}")

        i += 1

    return descCDK




def calculateDescriptorsObabel(smiles):
    descriptorsObabel = pd.DataFrame()
    columns = list(pybel.readstring("smi","C").calcdesc().keys())
    i = 0
    for smile in smiles:
        try:
            mol = pybel.readstring("smi",smile)
            desc = mol.calcdesc()
            descriptorsObabel = pd.concat([descriptorsObabel,pd.DataFrame(desc,index=[i])])
        except:
            descriptorsObabel = pd.concat([descriptorsObabel,pd.DataFrame(pd.NA,index=[i],columns=columns)])
            print(f"Obabel: Failed to calculate descriptors for {smi}")
        i+=1
    descriptorsObabel.drop(columns=["cansmi","cansmiNS","formula","title","InChI","InChIKey","smarts"],inplace=True)
    columns = [col + "_Obabel" for col in descriptorsObabel.columns]
    descriptorsObabel.columns = columns
    return descriptorsObabel

def calculateDescriptors(smiles):
    smilesCanon = []
    for smile in smiles:
        try:
            smilesCanon.append(CanonSmiles(smile))
        except:
            print(f"Failed to canonize {smile}")
            smilesCanon.append(smile)

    descriptorsCDK = calculateDescriptorsCDK(smilesCanon)
    descriptorsObabel = calculateDescriptorsObabel(smilesCanon)
    descriptorsRDKit = calculateDescriptorsRDKit(smilesCanon)
    descriptorsJazzy = calculateDescriptorsJazzy(smilesCanon)
    descriptorsConj = calculateConjugatedDescriptors(smilesCanon)
    descriptors = pd.concat([descriptorsJazzy,descriptorsCDK,descriptorsRDKit,descriptorsObabel,descriptorsConj],axis=1)
    descriptors.insert(0,"smiles",smilesCanon)
    return descriptors

def createDataset(data,fileName):
    data.sort_index(inplace=True)
    smiles_column = None
    if "smiles" in data.columns:
        smiles_column = "smiles"
    elif "SMILES" in data.columns:
        smiles_column = "SMILES"
    if smiles_column is not None:
        smiles = data[smiles_column].to_numpy()
        dataset = calculateDescriptors(smiles)
        for index in range(len(data.columns)-1,0,-1):
            dataset.insert(1,data.columns[index],data[data.columns[index]].to_numpy())
        dataset.to_csv("desc_"+fileName,index=False)

def main():
    reader = ParameterReader()
    parameters = reader.readParameters()
    fileName = parameters["data"]
    createDataset(pd.read_csv(fileName),fileName)

if __name__ == "__main__":
    main()
