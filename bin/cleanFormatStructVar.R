#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
struct_cnv_variants <- args[1]
struct_sv_variants <- args[2]
gene_name       <- args[3]


# Load annotation
annotation <- read.table(struct_cnv_variants, header = TRUE, sep = "\t")

# Do your cleaning / processing
cleaned_table <- annotation


# Output file
out_file <- paste0(gene_name, "_structural_variants.tsv")

# Write the clean table only
write.table(cleaned_table, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
