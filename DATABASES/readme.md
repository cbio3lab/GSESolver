# AqSolDBc.csv
Latest version of AQSolDB used for this work ([doi.org/10.1038/s41597-019-0151-1](https://www.nature.com/articles/s41597-019-0151-1)).

# AqSolDBc_noions.csv
AQSolDB after filtering every molecule with ionizable groups.

# Avdeef_sol_db.csv
Wiki-pS0 database extracted from *Predicting Solubility of New Drugs: Handbook of Critically Curated Data for Pharmaceutical Research* (**2024**).

# Avdeef_sol_db_298K.csv
Wiki-pS0 database after filtering only solubility values for each molecules closest to 298 K.

# BigSolDBv2.0.csv
Latest version of BigSolDB used for this work ([10.1038/s41597-025-05559-8](https://www.nature.com/articles/s41597-025-05559-8)).

# BigSolDBv2.0_w_298K.csv
BigSolDB2.0 after filtering only acqueous solubility values closest to 298 K.

# BigSolDBv2.0_w_298K_noions.csv
BigSolDB2.0 after filtering every molecule with ionizable groups.

# SolCbio3Database_merged.csv
SolCBio3Database after merging primary literature values with the external databases. Sources of each $\log{S_0}$, $mp$, and $\log{P_{\text{N}}}$ values are reported.

# descriptors_SolCbio3Database_merged.csv
Descriptors calculated for SolCbio3Database

# filtered_descriptors_SolCbio3Database_merged.csv
Descriptors for SolCbio3Database after UFS.

# RFE_SolCbio3_descs_fps.csv
Filtered descriptors for SolCbio3Database and fingerprints afer RFE.

# logPN_final.csv
In-house $\log{P_{\text{N}}}$ values database ([10.5281/zenodo.18616059](https://zenodo.org/records/18616059)).

# filter_NAs_db.R
Small script to filter empty values.

# ionizable_groups.py
Script used to filter molecules with ionizable groups

# merge_logs_datasets.py
Script used to merge datasets.
