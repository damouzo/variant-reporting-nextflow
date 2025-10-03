#!/usr/bin/env Rscript
# Script: cleanFormatBCIPatientsVar

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
gene_name        <- args[1]  # "DHX34"
vep_annot_file   <- args[2]  # "/mnt/c/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/results/DHX34/DHX34_bci_patients_vep_output.txt"
variant_file     <- args[3]  # 

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
  rename(Uploaded_variation = `#Uploaded_variation`) %>%      
  mutate(
    IMPACT = factor(IMPACT, ordered = TRUE,
            levels = c("HIGH", "MODERATE", "LOW", "MODIFIER")))

# Clean original data
variant_table <- variant_table %>%
    mutate(
        Uploaded_variation = ifelse(is.na(Transcript) | Transcript == "-",  NA,
                             paste0(Transcript, ":", cDNA_Change_HGVS_c)))


# Join vep annotations with patient data ---------------------------------------
annot_var_table <- variant_table %>%
  left_join(vep_table, by = "Uploaded_variation") %>%
  rename(Variant = Uploaded_variation) %>%
  distinct()

# Clean output table ---------------------------------------------------------------
annot_var_table <- annot_var_table %>%
  filter(CANONICAL == "YES") %>%
  filter(SYMBOL == gene_name)

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
