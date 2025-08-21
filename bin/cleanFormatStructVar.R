#!/usr/bin/env Rscript
# Script: cleanFormatStructVar

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
struct_cnv_variants <- args[1]
struct_sv_variants  <- args[2]
gene_name           <- args[3]

# Libraries  -------------------------------------------------------------------
library(tidyverse)

# Settings ----------------------------------------------------------------------
set.seed(23)


# Load CNV and SV files --------------------------------------------------------
CNV <- readRDS(struct_cnv_variants)
SV <- readRDS(struct_sv_variants)


# Combine Data -------------------------------------------------------------------
## Identify unique columns in each dataset
unique_SV_col <- colnames(SV)[!colnames(SV) %in% colnames(CNV) ] 
colnames(SV)[colnames(SV) %in% unique_SV_col] <- paste0(unique_SV_col, "_SVonly")

unique_CNV_col <- colnames(CNV)[!colnames(CNV) %in% colnames(SV) ] 
colnames(CNV)[colnames(CNV) %in% unique_CNV_col] <- paste0(unique_CNV_col, "_CNVonly")

CNV_SV <- bind_rows(CNV, SV)


# Outputs -----------------------------------------------------------------------
# Save rds file
output_file_rds <- paste0(gene_name, "_structural_variants.rds")
saveRDS(CNV_SV, output_file_rds)

# save tsv file
output_file_tsv <- paste0(gene_name, "_structural_variants.tsv")
write_tsv(CNV_SV, output_file_tsv)

# Print summary information
cat("Processing CNV SV combination for Gene:", gene_name, "\n")
cat("Output files:", output_file_rds, output_file_tsv, "\n")


