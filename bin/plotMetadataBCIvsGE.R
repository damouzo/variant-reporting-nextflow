#!/usr/bin/env Rscript
# Script: plotQCsmallVar

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: plotMetadataBCIvsGE.R <comparison_rds> <bci_clean_rds> <part_metadata_file> <gene_name>")
}

comparison_file     <- args[1]
bci_clean_file      <- args[2]
metadata_file       <- args[3]
gene_name           <- args[4]

# Message validation files
cat("Validating input files for gene:", gene_name, "\n")
missing_files <- c()
if (!file.exists(comparison_file)) missing_files <- c(missing_files, paste("Comparison file:", comparison_file))
if (!file.exists(bci_clean_file)) missing_files <- c(missing_files, paste("BCI clean file:", bci_clean_file))
if (!file.exists(metadata_file)) missing_files <- c(missing_files, paste("Metadata file:", metadata_file))

if (length(missing_files) > 0) {
  cat("ERROR: Missing files detected:\n")
  cat(paste(missing_files, collapse = "\n"), "\n")
  stop("Cannot proceed with missing input files")
}

cat("All input files validated successfully\n")

# Libraries  -------------------------------------------------------------------
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(paletteer)
library(grid)


# Settings ----------------------------------------------------------------------
set.seed(23)



###########
# Script a desarrollar
###########


# Create a simple overview table
