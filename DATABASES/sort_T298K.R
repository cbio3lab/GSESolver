setwd("C:/Users/ebert/OneDrive - Universidad de Costa Rica/SOLUBILIDAD/DATABASES")

# Enhanced version with detailed reporting
process_solubility_data_detailed <- function(data) {
  # Check if required columns exist
  required_cols <- c("SMILES", "logS0_t", "t")
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Calculate absolute difference from 25°C
  data$temp_diff <- abs(data$t - 25)
  
  # Count duplicates before processing
  smiles_count <- table(data$SMILES)
  duplicated_smiles <- names(smiles_count[smiles_count > 1])
  
  # Group by SMILES and keep only the row with minimum temperature difference
  result <- data %>%
    group_by(SMILES) %>%
    arrange(temp_diff) %>%
    slice(1) %>%
    ungroup() %>%
    select(-temp_diff)
  
  # Report what was done
  cat("=== Solubility Data Processing Report ===\n")
  cat("Total molecules (unique SMILES):", length(unique(data$SMILES)), "\n")
  cat("Molecules with multiple measurements:", length(duplicated_smiles), "\n")
  cat("Original number of measurements:", nrow(data), "\n")
  cat("Final number of measurements:", nrow(result), "\n")
  cat("Measurements removed:", nrow(data) - nrow(result), "\n")
  
  if (length(duplicated_smiles) > 0) {
    cat("\nMolecules with multiple measurements that were processed:\n")
    for (smiles in duplicated_smiles[1:min(5, length(duplicated_smiles))]) {
      orig_measurements <- data[data$SMILES == smiles, ]
      kept_measurement <- result[result$SMILES == smiles, ]
      cat(paste0("  ", smiles, ": ", nrow(orig_measurements), 
                 " measurements -> kept T = ", kept_measurement$t, "°C\n"))
    }
    if (length(duplicated_smiles) > 5) {
      cat("  ... and", length(duplicated_smiles) - 5, "more\n")
    }
  }
  
  return(result)
}


library(tidyverse)

data <- read.csv('Avdeef_sol_db.csv')
# Usage example:
processed_data <- process_solubility_data_detailed(data)


write_csv(processed_data, 'Avdeef_sol_db_298K.csv')
