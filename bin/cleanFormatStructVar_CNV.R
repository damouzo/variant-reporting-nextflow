#!/usr/bin/env Rscript

# Load required libraries
library(tidyverse)


# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
cnv_germline_file <- args[1]
cnv_somatic_file <- args[2]
gene_name <- args[3]

output_file <- paste0(gene_name, "_struct_cnv_variants.tsv")  # Generamos el nombre del archivo de salida

# Read and process germline CNV data
cnv_germline <- read_tsv(cnv_germline_file, show_col_types = FALSE) %>%
    mutate(origin = "germline")

# Read and process somatic CNV data
cnv_somatic <- read_tsv(cnv_somatic_file, show_col_types = FALSE) %>%
    mutate(origin = "somatic")

# Combine and process the data
combined_cnv <- bind_rows(cnv_germline, cnv_somatic) %>%
    # Add any additional processing steps here
    mutate(
        gene = gene_name,
        variant_type = "CNV"
    )

# Write the output
write_tsv(combined_cnv, output_file)
