#!/usr/bin/env Rscript
# Script: cleanFormatStructVar_SV

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
sv_germline_file <- args[1] # "c:/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/data/raw_data/WGS_Variants/structuralVar/SV_chr19_DHX34_ENSG223423_GRCh38_germline.tsv"
sv_somatic_file <- args[2] # "c:/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/data/raw_data/WGS_Variants/structuralVar/SV_chr19_DHX34_ENSG234234_GRCh38_somatic.tsv"
gene_name <- args[3] # "DHX34"

# Libraries  -------------------------------------------------------------------
library(tidyverse)


# Settings ----------------------------------------------------------------------
set.seed(23)

# Load SV files -----------------------------------------------------------------
sv_germline <- read_tsv(sv_germline_file, show_col_types = FALSE)
sv_somatic <- read_tsv(sv_somatic_file, show_col_types = FALSE)



# Clean Data -------------------------------------------------------------------
## Add variant Origin
sv_germline <- sv_germline %>% mutate(variant_origin = "germline")
sv_somatic <- sv_somatic %>% mutate(variant_origin = "somatic")

# Solve GT and GQ diff
if(!"GT" %in% colnames(sv_somatic)) {sv_somatic$GT <- NA}
if(!"GQ" %in% colnames(sv_somatic)) {sv_somatic$GQ <- NA}

# Combine clean columns
SV <- rbind(sv_germline, sv_somatic)
colnames(SV)[1] <- c("SAMPLE")

# Format Columns
numeric_columns <- c("QUAL","INFO_END", "OVERLAP_CNT", "OVERLAP_PCT", "INFO_SVLEN")
SV <- SV %>% mutate(across(any_of(numeric_columns), as.numeric))

# Generate columns
SV$algorithm <- str_sub(sapply(strsplit(SV$ID, split=":"), `[`, 1), end = 5)
SV$SYMBOL_annotation <- gene_name
SV$variant_origin_file <- "SV"


# Output -----------------------------------------------------------------------
# Write the output file
output_file <- paste0(gene_name, "_struct_sv_variants.rds")
saveRDS(SV, output_file)

# Print summary information
cat("Processing SVs completed successfully for Gene:", gene_name, "\n")
cat("Germline variants:", sum(SV$variant_origin == "germline", na.rm = TRUE), "\n")
cat("Somatic variants:", sum(SV$variant_origin == "somatic", na.rm = TRUE), "\n")
cat("Output file:", output_file, "\n")


