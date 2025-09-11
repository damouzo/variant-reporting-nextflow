#!/usr/bin/env Rscript
# Script: cleanFormatBCIPatientsVar

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
vep_annot_file   <- args[1]  # "/mnt/c/scratch/2e/56b66f4d4935c698b3dac8bc29f26b/DDX41_vep_output.txt"
gene_name        <- args[2]  # "DDX41"
variant_file     <- args[3]  # "/mnt/c/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/data/raw_data/BCI_patients/250829_DDX41_patients_var_data.tsv"

# Libraries  -------------------------------------------------------------------
library(tidyverse)


# Settings ----------------------------------------------------------------------
set.seed(23)


# Load gene vep annotated and original variant file ------------------------------
# Load original data
variant_table <- read_tsv(variant_file)

# Load VEP data
vep_table <- read_tsv(vep_annot_file, comment = "##", na = c("-", ""),
                      col_types = cols(Location = col_character(),
                                        .default = col_guess()))


# Clean tables -------------------------------------------------------------------
# Vep table output
vep_table <- vep_table %>%
  rename(Uploaded_variation = `#Uploaded_variation`) %>%        # quita el '#'
  mutate(
    IMPACT = factor(IMPACT, ordered = TRUE,
            levels = c("HIGH", "MODERATE", "LOW", "MODIFIER")))

# Clean original data
variant_table <- variant_table %>%
    mutate(
        Uploaded_variation = ifelse(is.na(Transcript) | Transcript == "-",  NA,
                             paste0(gsub("\\..*", ":", Transcript), cDNA_Change_HGVS_c.)))


# Join vep annotations with patient data ---------------------------------------
annot_var_table <- vep_table %>%
  inner_join(variant_table, by = "Uploaded_variation") %>%
  rename(Variant = Uploaded_variation)


# Print summary ----------------------------------------------------------------
cat("\nGlimpse of table data of Variants \n")
glimpse(annot_var_table)


# Save variant table  --------------------------------------------------------------
# TSV
out_file_tsv <- paste0(gene_name, "_bci_patients_variants.tsv")
write.table(annot_var_table, file = out_file_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# RDS
out_file_rds <- paste0(gene_name, "_bci_patients_variants.rds")
saveRDS(annot_var_table, file = out_file_rds)

# Print output messages
cat("\nOutput files written:")
cat("\n- TSV file:", out_file_tsv)
cat("\n- RDS file:", out_file_rds, "\n")
