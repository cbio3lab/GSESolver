from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors
import pandas as pd

from rdkit.Chem import GraphDescriptors
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski



def calculateDescriptorsMolRDKitOriginal(molecule):
    descriptors = {}

    descriptors["ExactMolWt"] = Descriptors.ExactMolWt(molecule)
    descriptors["FpDensityMorgan1"] = Descriptors.FpDensityMorgan1(molecule)
    descriptors["FpDensityMorgan2"] = Descriptors.FpDensityMorgan2(molecule)
    descriptors["FpDensityMorgan3"] = Descriptors.FpDensityMorgan3(molecule)

    descriptors["HeavyAtomMolWt"] = Descriptors.HeavyAtomMolWt(molecule)

    descriptors["MaxAbsPartialCharge"] = Descriptors.MaxAbsPartialCharge(molecule)
    descriptors["MaxPartialCharge"] = Descriptors.MaxPartialCharge(molecule)
    descriptors["MinAbsPartialCharge"] = Descriptors.MinAbsPartialCharge(molecule)
    descriptors["MinPartialCharge"] = Descriptors.MinPartialCharge(molecule)

    descriptors["MolWt"] = Descriptors.MolWt(molecule)
    descriptors["NumRadicalElectrons"] = Descriptors.NumRadicalElectrons(molecule)
    descriptors["NumValenceElectrons"] = Descriptors.NumValenceElectrons(molecule)

    descriptors["BalabanJ"] = GraphDescriptors.BalabanJ(molecule)
    descriptors["BertzCT"] = GraphDescriptors.BertzCT(molecule)

    descriptors["Chi0"] = GraphDescriptors.Chi0(molecule)
    descriptors["Chi1"] = GraphDescriptors.Chi1(molecule)

    descriptors["Chi0n"] = GraphDescriptors.Chi0n(molecule)
    descriptors["Chi1n"] = GraphDescriptors.Chi1n(molecule)
    descriptors["Chi2n"] = GraphDescriptors.Chi2n(molecule)
    descriptors["Chi3n"] = GraphDescriptors.Chi3n(molecule)
    descriptors["Chi4n"] = GraphDescriptors.Chi4n(molecule)

    descriptors["Chi0v"] = GraphDescriptors.Chi0v(molecule)
    descriptors["Chi1v"] = GraphDescriptors.Chi1v(molecule)
    descriptors["Chi2v"] = GraphDescriptors.Chi2v(molecule)
    descriptors["Chi3v"] = GraphDescriptors.Chi3v(molecule)
    descriptors["Chi4v"] = GraphDescriptors.Chi4v(molecule)

    descriptors["HallKierAlpha"] = GraphDescriptors.HallKierAlpha(molecule)

    descriptors["Kappa1"] = GraphDescriptors.Kappa1(molecule)
    descriptors["Kappa2"] = GraphDescriptors.Kappa2(molecule)
    descriptors["Kappa3"] = GraphDescriptors.Kappa3(molecule)

    descriptors["MolLogP"] = Crippen.MolLogP(molecule)
    descriptors["MolMR"] = Crippen.MolMR(molecule)

    descriptors["HeavyAtomCount"] = Lipinski.HeavyAtomCount(molecule)
    descriptors["NHOHCount"] = Lipinski.NHOHCount(molecule)
    descriptors["NOCount"] = Lipinski.NOCount(molecule)
    descriptors["NumHAcceptors"] = Lipinski.NumHAcceptors(molecule)
    descriptors["NumHDonors"] = Lipinski.NumHDonors(molecule)
    descriptors["NumHeteroatoms"] = Lipinski.NumHeteroatoms(molecule)
    descriptors["NumRotatableBonds"] = Lipinski.NumRotatableBonds(molecule)
    descriptors["NumAmideBonds"] = rdMolDescriptors.CalcNumAmideBonds(molecule)

    descriptors["NumAromaticRings"] = Lipinski.NumAromaticRings(molecule)
    descriptors["NumAliphaticRings"] = Lipinski.NumAliphaticRings(molecule)
    descriptors["NumSaturatedRings"] = Lipinski.NumSaturatedRings(molecule)

    descriptors["NumAliphaticCarbocycles"] = Lipinski.NumAliphaticCarbocycles(molecule)
    descriptors["NumAliphaticHeterocycles"] = Lipinski.NumAliphaticHeterocycles(molecule)

    descriptors["NumAromaticCarbocycles"] = Lipinski.NumAromaticCarbocycles(molecule)
    descriptors["NumAromaticHeterocycles"] = Lipinski.NumAromaticHeterocycles(molecule)

    descriptors["NumSaturatedCarbocycles"] = Lipinski.NumSaturatedCarbocycles(molecule)
    descriptors["NumSaturatedHeterocycles"] = Lipinski.NumSaturatedHeterocycles(molecule)

    descriptors["RingCount"] = Lipinski.RingCount(molecule)
    descriptors["FractionCSP3"] = Lipinski.FractionCSP3(molecule)

    descriptors["NumSpiroAtoms"] = rdMolDescriptors.CalcNumSpiroAtoms(molecule)
    descriptors["NumBridgeheadAtoms"] = rdMolDescriptors.CalcNumBridgeheadAtoms(molecule)
    descriptors["TPSA"] = rdMolDescriptors.CalcTPSA(molecule)
    descriptors["LabuteASA"] = rdMolDescriptors.CalcLabuteASA(molecule)
    
    PEOE_VSA = rdMolDescriptors.PEOE_VSA_(molecule)
    for i in range(1,15):
        descriptors["PEOE_VSA"+str(i)] = PEOE_VSA[i-1]

    SMR_VSA = rdMolDescriptors.SMR_VSA_(molecule)
    for i in range(1,11):
        descriptors["SMR_VSA"+str(i)] = SMR_VSA[i-1]

    SlogP_VSA = rdMolDescriptors.SlogP_VSA_(molecule)
    for i in range(1,11):
        descriptors["SlogP_VSA"+str(i)] = SlogP_VSA[i-1]
    
    
    MQNs = rdMolDescriptors.MQNs_(molecule)
    for i in range(1,43):
        descriptors["MQNs_"+str(i)] = MQNs[i-1]

    return descriptors

def calculateDescriptorsMolRDKit(molecule):
    descriptors = Descriptors.CalcMolDescriptors(molecule)
    descriptors["NumAmideBonds"] = rdMolDescriptors.CalcNumAmideBonds(molecule)
    descriptors["NumSpiroAtoms"] = rdMolDescriptors.CalcNumSpiroAtoms(molecule)
    descriptors["NumBridgeheadAtoms"] = rdMolDescriptors.CalcNumBridgeheadAtoms(molecule)
    PEOE_VSA = rdMolDescriptors.PEOE_VSA_(molecule)
    for i in range(1,15):
        descriptors["PEOE_VSA_"+str(i)] = PEOE_VSA[i-1]

    SMR_VSA = rdMolDescriptors.SMR_VSA_(molecule)
    for i in range(1,11):
        descriptors["SMR_VSA_"+str(i)] = SMR_VSA[i-1]

    SlogP_VSA = rdMolDescriptors.SlogP_VSA_(molecule)
    for i in range(1,11):
        descriptors["SlogP_VSA_"+str(i)] = SlogP_VSA[i-1]
    
    
    MQNs = rdMolDescriptors.MQNs_(molecule)
    for i in range(1,43):
        descriptors["MQNs_"+str(i)] = MQNs[i-1]
    return descriptors
    


def calculateDescriptorsRDKit(smiles):
    descriptorsRDKit = {}
    for smile in smiles:
        try:
            molecule = Chem.MolFromSmiles(smile)
            descriptors = calculateDescriptorsMolRDKit(molecule)
            for descriptor in descriptors:
                if (descriptor+"_rdkit") in descriptorsRDKit:
                    descriptorsRDKit[descriptor+"_rdkit"].append(descriptors[descriptor])
                else:
                    descriptorsRDKit[descriptor+"_rdkit"] = [descriptors[descriptor]]
        except:
            print(f"RDKit: Faile to calculate descriptors for {smile}")
            for descriptor in descriptorsRDKit:
                descriptorsRDKit[descriptor].append("")
    return pd.DataFrame(descriptorsRDKit)
