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


# Filter for canonical variants -------------------------------------------------
# Stats of variant filtering
stats_variants_filter <- data.frame(metric = character(), count = numeric(), stringsAsFactors = FALSE)
stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_raw_variants", 
                                        count = nrow(variant_table)))

# Filter for canonical variants
variant_filtered_table <- variant_table[variant_table$CANONICAL_annotation == "YES", ]

# Add canonical variants count to stats
stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_canonical_variants", 
                                        count = nrow(variant_filtered_table)))


# Filter ClinVar annotations ---------------------------------------------------
# Function to categorize ClinVar annotations
categorize_clinvar <- function(annotation) {
  if (is.na(annotation) || annotation == "-") { return("-") }
  
  # Split by comma and slash, then trim whitespace
  terms <- unlist(strsplit(annotation, "[,/]"))
  terms <- trimws(terms)
  
  # Define groups
  pathogenic_terms <- c("pathogenic", "likely_pathogenic")
  vus_terms <- c("uncertain_significance", "conflicting_interpretations_of_pathogenicity")
  benign_terms <- c("benign", "likely_benign")
  if (any(terms %in% pathogenic_terms)) {return("pathogenic") }
  if (any(terms %in% vus_terms)) {return("VUS")}
  if (any(terms %in% benign_terms)) {return("benign")}
  return("VUS") # Default to VUS for other cases
}

# Apply categorization
variant_filtered_table$clinvar_category <- sapply(variant_filtered_table$CLIN_SIG_annotation, categorize_clinvar)

# Count each category and add to stats
clinvar_counts <- table(variant_filtered_table$clinvar_category)
for (category in names(clinvar_counts)) {
  stats_variants_filter <- rbind(stats_variants_filter, 
                                data.frame(metric = paste0("total_canonical_", category), 
                                          count = as.numeric(clinvar_counts[category])))
}

# Filter out benign variants
variant_filtered_table <- variant_filtered_table[variant_filtered_table$clinvar_category != "benign", ]

# Add count after benign filtering
stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_canonical_nonBenign", 
                                        count = nrow(variant_filtered_table)))


# Filter by MAX_AF annotation --------------------------------------------------
# Keep variants with MAX_AF = "-" or < 0.001
variant_filtered_table <- variant_filtered_table[
  variant_filtered_table$MAX_AF_annotation == "-" | 
  (variant_filtered_table$MAX_AF_annotation != "-" & 
   as.numeric(variant_filtered_table$MAX_AF_annotation) < 0.001), ]

# Add count after MAX_AF filtering
stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_canonical_nonBenign_MAF<0.001", 
                                        count = nrow(variant_filtered_table)))


# Print summary ----------------------------------------------------------------
cat("\nGlimpse of table data of Variants \n")
glimpse(variant_filtered_table)


# Save variant table  --------------------------------------------------------------
# Save stats file
stats_file_csv <- paste0(gene_name, "_small_variants_filtered_stats.csv")
write.csv(stats_variants_filter, file = stats_file_csv, row.names = FALSE)

# TSV
out_file_tsv <- paste0(gene_name, "_small_variants_filtered.tsv")
write.table(variant_filtered_table, file = out_file_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# RDS
out_file_rds <- paste0(gene_name, "_small_variants_filtered.rds")
saveRDS(variant_filtered_table, file = out_file_rds)

# Print output messages
cat("\nOutput files written:")
cat("\n- Stats file:", stats_file_csv)
cat("\n- TSV file:", out_file_tsv)
cat("\n- RDS file:", out_file_rds, "\n")
