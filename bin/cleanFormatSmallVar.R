#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
annotation_file <- args[1]
gene_name       <- args[2]


# Load annotation
annotation <- read.table(annotation_file, header = TRUE, sep = "\t")

# Do your cleaning / processing
cleaned_table <- annotation


# Output file
out_file <- paste0(gene_name, "_small_variants.tsv")

# Write the clean table only
write.table(cleaned_table, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
