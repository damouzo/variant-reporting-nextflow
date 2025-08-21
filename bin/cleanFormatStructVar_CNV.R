#!/usr/bin/env Rscript
# Script: cleanFormatStructVar_CNV

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
cnv_germline_file <- args[1] # "C:/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/data/raw_data/WGS_Variants/structuralVar/CNV_chr5_DDX41_ENSG234234_GRCh38_germline.tsv"
cnv_somatic_file <- args[2] # "C:/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/data/raw_data/WGS_Variants/structuralVar/CNV_chr5_DDX41_ENSG234234_GRCh38_somatic.tsv"
gene_name <- args[3] # "DDX41"

# Libraries  -------------------------------------------------------------------
library(tidyverse)

# Settings ----------------------------------------------------------------------
set.seed(23)

# Load SV files -----------------------------------------------------------------
cnv_germline <- read_tsv(cnv_germline_file, show_col_types = FALSE)
cnv_somatic <- read_tsv(cnv_somatic_file, show_col_types = FALSE)


# Clean Data -------------------------------------------------------------------
## Add variant Origin 
cnv_germline <- cnv_germline %>% mutate(variant_origin = "germline")
cnv_somatic <- cnv_somatic %>% mutate(variant_origin = "somatic")

# Solve GT diff
if(!"GT" %in% colnames(cnv_somatic)) {cnv_somatic$GT <- NA}

# Combine clean columns
CNV <- rbind(cnv_germline, cnv_somatic)
colnames(CNV)[1] <- c("SAMPLE")

# Format Columns
numeric_columns <- c("MCC")
CNV <- CNV %>% mutate(across(all_of(numeric_columns), as.numeric))

# Generate columns
CNV$algorithm <- sapply(strsplit(CNV$ID, split=":"), `[`, 1)
CNV$SYMBOL_annotation <- gene_name
CNV$variant_origin_file <- "CNV"



# Output -----------------------------------------------------------------------
# Write the output file
output_file <- paste0(gene_name, "_struct_cnv_variants.rds")
saveRDS(CNV, output_file)

# Print summary information
cat("Processing CNVs completed successfully for Gene:", gene_name, "\n")
cat("Germline variants:", sum(CNV$variant_origin == "germline", na.rm = TRUE), "\n")
cat("Somatic variants:", sum(CNV$variant_origin == "somatic", na.rm = TRUE), "\n")
cat("Output file:", output_file, "\n")
