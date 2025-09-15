#!/usr/bin/env Rscript
# Script: compareBCIwithGEvar

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
bci_file <- args[1]  # BCI small variants RDS file
ge_file <- args[2]    # GE small variants RDS file
gene_name <- args[3]  # Gene name


# Libraries  -------------------------------------------------------------------
library(tidyverse)


# Settings ----------------------------------------------------------------------
set.seed(23)


# Load gene annotated variant file ---------------------------------------------
bci_var_table <- readRDS(bci_file)
ge_variant_table <- readRDS(ge_file)


##############
# Script to compare BCI patients variants with GE variants for a specific gene
##############

# Print summary ----------------------------------------------------------------
cat("\nGlimpse of table data of Variants \n")
glimpse(compared_table)


# Save variant table  --------------------------------------------------------------
# TSV
out_file_tsv <- paste0(gene_name, "compared_BCI_GE_small_variants.tsv")
write.table(compared_table, file = out_file_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# RDS
out_file_rds <- paste0(gene_name, "compared_BCI_GE_small_variants.rds")
saveRDS(compared_table, file = out_file_rds)

# Print output messages
cat("\nOutput files written:")
cat("\n- TSV file:", out_file_tsv)
cat("\n- RDS file:", out_file_rds, "\n")
