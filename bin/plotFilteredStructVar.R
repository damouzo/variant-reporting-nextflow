#!/usr/bin/env Rscript
# Script: plotFilteredStructVar

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("script.R <filtered_table> <gene_name> <prot_file> <exon_file> <part_metadata_file>")
}

filtered_table  <- args[1] # "C:\\Users\\qp241615\\OneDrive - Queen Mary, University of London\\Documents\\4. Projects\\1. DHX34\\results\\DHX34\\DHX34_small_variants.rds"
gene_name       <- args[2] # "DHX34"
prot_file       <- args[3] # "C:\\Users\\qp241615\\OneDrive - Queen Mary, University of London\\Documents\\4. Projects\\1. DHX34\\data\\reference\\Protein\\DHX34.gff"
exon_file       <- args[4] # "C:\\Users\\qp241615\\OneDrive - Queen Mary, University of London\\Documents\\4. Projects\\1. DHX34\\data\\reference\\Exon\\DHX34.tsv"
p_metadata_file <- args[5] # "C:\\Users\\qp241615\\OneDrive - Queen Mary, University of London\\Documents\\4. Projects\\1. DHX34\\results\\DHX34\\DHX34_small_variants_participantMetadata.tsv"

# Message validation files
cat("Validating input files for gene:", gene_name, "\n")
missing_files <- c()
if (!file.exists(filtered_table)) missing_files <- c(missing_files, paste("Filtered table:", filtered_table))
if (!file.exists(prot_file)) missing_files <- c(missing_files, paste("Protein file:", prot_file))
if (!file.exists(exon_file)) missing_files <- c(missing_files, paste("Exon file:", exon_file))
if (!file.exists(p_metadata_file)) missing_files <- c(missing_files, paste("Participant metadata file:", p_metadata_file))

if (length(missing_files) > 0) {
  cat("ERROR: Missing files detected:\n")
  cat(paste(missing_files, collapse = "\n"), "\n")
  stop("Cannot proceed with missing input files")
}

cat("All input files validated successfully\n")


# Libraries  -------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(circlize)
library(paletteer)


# Settings ----------------------------------------------------------------------
set.seed(23)


# Load Data ---------------------------------------------------------------------
## Structural Variants
variants_table <- readRDS(struct_variants)

## Participant Metadata
metadata_info <- readRDS(part_metadata_file)





cat("Plots generated successfully for gene:", gene_name, "\n")