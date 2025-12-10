library(tidyverse)

sol <- read.csv('SolCbio3Database_merged.csv')

sol_sinlogPNA <- sol %>% filter(!is.na(logPN))
sol_sinmpNA <- sol_sinlogPNA %>% filter(!is.na(mp))

#write.csv(sol_sinmpNA, file = "SolCbio3Database_merged_noNA.csv")



# Remove duplicate molecules based on SMILES codes
remove_duplicate_molecules <- function(data) {
  # Check if required column exists
  if (!"smiles" %in% colnames(data)) {
    stop("Required column 'smiles' not found in the dataframe")
  }
  
  # Count duplicates before processing
  original_rows <- nrow(data)
  unique_smiles_before <- length(unique(data$smiles))
  duplicate_count <- sum(duplicated(data$smiles))
  
  # Remove duplicate rows based on SMILES column
  # keep = "first" keeps the first occurrence of each duplicate
  cleaned_data <- data[!duplicated(data$smiles), ]
  
  # Alternative using dplyr (uncomment if you prefer this approach)
  # library(dplyr)
  # cleaned_data <- data %>% distinct(smiles, .keep_all = TRUE)
  
  # Report results
  cat("=== Duplicate Removal Report ===\n")
  cat("Original number of rows:", original_rows, "\n")
  cat("Final number of rows:", nrow(cleaned_data), "\n")
  cat("Rows removed:", original_rows - nrow(cleaned_data), "\n")
  cat("Unique molecules before:", unique_smiles_before, "\n")
  cat("Unique molecules after:", length(unique(cleaned_data$smiles)), "\n")
  cat("Duplicate SMILES found:", duplicate_count, "\n")
  
  return(cleaned_data)
}


sol_unique <- remove_duplicate_molecules(sol)
write.csv(sol_sinmpNA, file = "SolCbio3Database_merged_noNA.csv")
