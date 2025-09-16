#!/usr/bin/env Rscript
# Script: filterSmallVar

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
clean_smallvar_file <- args[1]  # "C:/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/data/raw_data/WGS_Variants/smallVar/GRCh38_DDX41_annotated_variants.tsv"
gene_name <- args[2]        # "DDX41"

# Libraries  -------------------------------------------------------------------
library(tidyverse)


# Settings ----------------------------------------------------------------------
set.seed(23)


# Load gene annotated variant file ---------------------------------------------
variant_table <- readRDS(clean_smallvar_file)

#####################
# Script to filter variants based on specific criteria
# Load critearia from a config file or define them here? 
# MAF -> MAX_AF ≤ 0.01 (1%)
####################

# Print summary ----------------------------------------------------------------
cat("\nGlimpse of table data of Variants \n")
glimpse(variant_table)


# Save variant table  --------------------------------------------------------------
# TSV
out_file_tsv <- paste0(gene_name, "_small_variants_filtered.tsv")
write.table(variant_table, file = out_file_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# RDS
out_file_rds <- paste0(gene_name, "_small_variants_filtered.rds")
saveRDS(variant_table, file = out_file_rds)

# Print output messages
cat("\nOutput files written:")
cat("\n- TSV file:", out_file_tsv)
cat("\n- RDS file:", out_file_rds, "\n")
