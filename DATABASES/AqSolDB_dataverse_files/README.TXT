File AqSolDBc.csv
Curated data from the AqSolDB. The available columns are:
-* ID: Compound ID (string)
-* InChI: InChI code of the chemical structure (string)
-* Solubility Mole/L: logarithm of the thermodynamic solubility in water at pH 7 (+/-1) at ~300K (float)
-* SMILEScurated: Curated SMILES code of the chemical structure (string)
-* SD: Standard laboratory Deviation, default value: -1 (float)
-* Group: Data quality label imported from AqSolDB (string)
-* Dataset: Source of the data point (string)
-* Composition: Purity of the substance: mono-constituent, multi-constituent, UVCB (Categorical)
-* Error: Identifier error on the data point, default value: None (String)
-* Charge: Estimated formal charge of the compound at pH 7: Positive, Negative, Zwiterion, Uncharged (Categorical)

File OChemUnseen.csv
Solubility data from OChem, curated and orthogonal to AqSolDB. The available columns are:
-* SMILES: Curated SMILES code of the chemical structure (string)
-* LogS Mole/L: logarithm of the thermodynamic solubility in water at pH 7 (+/-1) (float)

File OChemOverlapping.csv
Solubility data from OChem, curated; chemical structures are also present inside AqSolDB. The available columns are:
-* SMILES: Curated SMILES code of the chemical structure (string)
-* LogS Mole/L: logarithm of the thermodynamic solubility in water at pH 7 (+/-1) (float)

File OChemCurated.csv
Solubility data from OChem, curated. The available columns are:
-* ID: Compound ID (string)
-* Name: Compound name (string)
-* SMILES: Curated SMILES code of the chemical structure (string)
-* SDi: Standard laboratory Deviation, default value: -1 (float)
-* Reference: Unformated bibliographic reference which the data point is originating from (string)
-* LogS Mole/L: logarithm of the thermodynamic solubility in water at pH 7 (+/-1) (float)
-* EXTERNALID: Compound ID as appearing in its data source, default value: None (string)
-* CASRN: CAS number of the compound, default value: None (string)
-* ARTICLEID: Source ID linked to the column Reference (string)
Temperature Temperature of the measure, in K (float)

