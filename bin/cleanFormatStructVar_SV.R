#!/usr/bin/env Rscript

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
sv_germline_file <- args[1]
sv_somatic_file <- args[2]
gene_name <- args[3]

# Define output file name
output_file <- paste0(gene_name, "_struct_sv_variants.tsv")

# Load required libraries
suppressPackageStartupMessages({
    library(tidyverse)
})

# Read and process germline SV data
sv_germline <- read_tsv(sv_germline_file, show_col_types = FALSE) %>%
    mutate(origin = "germline")

# Read and process somatic SV data
sv_somatic <- read_tsv(sv_somatic_file, show_col_types = FALSE) %>%
    mutate(origin = "somatic")

# Combine and process the data
combined_sv <- bind_rows(sv_germline, sv_somatic) %>%
    # Add any additional processing steps here
    mutate(
        gene = gene_name,
        variant_type = "SV"
    )

# Write the output
write_tsv(combined_sv, output_file)
