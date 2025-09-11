#!/usr/bin/env Rscript
# Script: cleanFormatBCIPatientsVar

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
variant_file  <- args[1]  # "C:/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/data/raw_data/BCI_patients/250829_DDX41_patients_var_data.tsv"
gene_name     <- args[2]  # "DDX41"

# Libraries  -------------------------------------------------------------------
library(tidyverse)


# Settings ----------------------------------------------------------------------
set.seed(23)


# Load gene annotated variant file ---------------------------------------------
variant_table <- read_tsv(variant_file)

#############################
# Script
############################

# Print summary ----------------------------------------------------------------
cat("\nGlimpse of table data of Variants \n")
glimpse(variant_table)


# Save variant table  --------------------------------------------------------------
# TSV
out_file_tsv <- paste0(gene_name, "_bci_patients_variants.tsv")
write.table(variant_table, file = out_file_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# RDS
out_file_rds <- paste0(gene_name, "_bci_patients_variants.rds")
saveRDS(variant_table, file = out_file_rds)

# Print output messages
cat("\nOutput files written:")
cat("\n- TSV file:", out_file_tsv)
cat("\n- RDS file:", out_file_rds, "\n")
