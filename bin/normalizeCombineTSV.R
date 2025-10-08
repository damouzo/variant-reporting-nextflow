#!/usr/bin/env Rscript
# Script: normalizeCombineTSV.R 

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript normalizeCombineTSV.R <gene_name> <file1.tsv> [file2.tsv] ...")
}

gene_name <- args[1]  # "DHX34"
input_files <- args[-1]  # c("/mnt/c/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/data/raw_data/BCI_patients/250926_DHX34_bci_patients_var_data.tsv", "/mnt/c/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/data/raw_data/BCI_patients/250926_DHX34_pv_patients_var_data.tsv")

# Libraries  -------------------------------------------------------------------
library(tidyverse)

# Settings ----------------------------------------------------------------------
set.seed(23)

# Functions --------------------------------------------------------------------
# Function to detect center from filename
detect_center <- function(filename) {
  filename_lower <- str_to_lower(basename(filename))
  if (str_detect(filename_lower, "_pv_")) return("PV")
  if (str_detect(filename_lower, "_bci_")) return("BCI") 
  if (str_detect(filename_lower, "_kc_")) return("KC")
  
  # If no pattern matches, try to infer from content later
  cat("DEBUG: No center pattern matched, returning UNKNOWN\n")
  return("UNKNOWN")
}

# Function to standardize column names based on center-specific patterns
normalize_columns <- function(df, center, filename) {
    # Create a copy to work with
    normalized_df <- df
    column_selected <- c("Variant_REF", "Patient_Centre_REF", "Gene_id", "Transcript", "cDNA_Change_HGVS_c", "Protein_Change_HGVS_p", "Center")
    
    # Center-specific column mapping based on actual data structures
    if (center == "BCI") {
        cat("\nApplying BCI center column mapping...\n")
        normalized_df <- df %>%
            mutate(
                Patient_Centre_REF = `FAMILY`,
                Gene_id = `GENE`,
                Transcript = `TRANSCRIPT`,
                cDNA_Change_HGVS_c = `VARIANT`,
                Protein_Change_HGVS_p = str_extract(cDNA_Change_HGVS_c, "p\\.\\([^)]+\\)"),
                cDNA_Change_HGVS_c = str_extract(cDNA_Change_HGVS_c, "c\\.[^:]+"),
                Center = "BCI"
            ) %>%
            mutate(Variant_REF = paste0("B", row_number())) %>%
            select(all_of(column_selected))
        
    } else if (center == "PV") {
        normalized_df <- df %>%
            mutate(
                Patient_Centre_REF = ID,
                Gene_id = GENE,
                Transcript = NM,
                cDNA_Change_HGVS_c = HGVS_C,
                Protein_Change_HGVS_p = HGVS_P,
                Reference_Allele = str_extract(cDNA_Change_HGVS_c, "(?<=\\d)[A-Z](?=>)"),
                Alternative_Allele = str_extract(cDNA_Change_HGVS_c, "(?<=>)[A-Z]"),
                Center = "PV"
            ) %>%
            mutate(Variant_REF = paste0("P", row_number())) %>%
            select(all_of(column_selected))
        
    } else if (center == "KC") {
        cat("\nApplying KC center column mapping...\n")
        normalized_df <- df %>%
            mutate(
                Patient_Centre_REF = Patient_REF,
                Gene_id = GenMut,
                Transcript = Transcript,
                cDNA_Change_HGVS_c = cDNA_Change_HGVS_c.,
                Protein_Change_HGVS_p = Protein_Change_HGVS_p.,
                Transcript = str_remove(Transcript, "\\..*$"),
                Reference_Allele = str_extract(cDNA_Change_HGVS_c, "(?<=\\d)[A-Z](?=>)"),
                Alternative_Allele = str_extract(cDNA_Change_HGVS_c, "(?<=>)[A-Z]"),
                Center = "KC"
            ) %>%
            mutate(Variant_REF = paste0("K", row_number())) %>%
            select(all_of(column_selected))
    } 
    return(normalized_df)
}


# Normalization processing --------------------------------------------------------------
combined_data <- tibble()

# Process each file
for (i in seq_along(input_files)) {
    # Read file
    df <- read_tsv(input_files[i])
    center <- detect_center(input_files[i])

    # Normalize columns
    normalized_df <- normalize_columns(df, center, input_files[i])
    
    # Add source file information
    normalized_df <- normalized_df %>%
        mutate(Source_File = basename(input_files[i]))
    
    # Combine with previous data
    combined_data <- bind_rows(combined_data, normalized_df)
}

# Save results ----------------------------------------------------------------

output_file <- paste0(gene_name, "_normalized_combined_data.tsv")
write_tsv(combined_data, output_file)

# Create VEP input file -------------------------------------------------------
vep_input_data <- combined_data %>%
    filter(!is.na(cDNA_Change_HGVS_c) & !is.na(Transcript)) %>%
    filter(cDNA_Change_HGVS_c != "NA" & Transcript != "NA") %>%
    mutate(
        variant_id = paste(Variant_REF, row_number(), sep="_"), # Create unique identifier
        transcript_clean = str_remove(Transcript, "\\..*"), # Clean transcript (remove version if exists)
        vep_input_line = paste0(transcript_clean, ":", cDNA_Change_HGVS_c) # Format for VEP: transcript:hgvs_c
    ) %>%
    select(variant_id, vep_input_line)

# Write VEP input file (without headers)
vep_output_file <- paste0(gene_name, "_vep_input.txt")
write_tsv(vep_input_data[,2], vep_output_file, col_names = FALSE)

cat("Prepared", nrow(vep_input_data), "variants from combined data for VEP annotation\n")

