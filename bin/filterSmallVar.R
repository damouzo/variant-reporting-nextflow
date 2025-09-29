#!/usr/bin/env Rscript
# Script: filterSmallVar

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
clean_smallvar_file   <- args[1]  # "/mnt/c/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/results/DHX34/DHX34_small_variants.rds
gene_name             <- args[2]  # "DHX34"
part_metadata_file    <- args[3]  # "/mnt/c/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/results/DHX34/DHX34_structural_variants_participantMetadata.rds"

# Libraries  -------------------------------------------------------------------
library(tidyverse)

# Settings ----------------------------------------------------------------------
set.seed(23)

# Load gene annotated variant file ---------------------------------------------
variant_table <- readRDS(clean_smallvar_file)
participant_metadata <- readRDS(part_metadata_file)


# Function to count unique participants ----------------------------------------
count_unique_participants <- function(variant_df, metadata_df) {
  # Extract all samples from Het_samples, Hom_samples, and Hemi_samples
  all_samples <- c()
  
  for (col in c("Het_samples", "Hom_samples", "Hemi_samples")) {
    if (col %in% colnames(variant_df)) {
      samples_in_col <- unlist(strsplit(variant_df[[col]], ","))
      samples_in_col <- trimws(samples_in_col)
      samples_in_col <- samples_in_col[samples_in_col != "NA" & samples_in_col != ""]
      all_samples <- c(all_samples, samples_in_col)
    }
  }
  
  # Get unique participant_ids
  unique_samples <- unique(all_samples)
  matched_participants <- metadata_df$participant_id[metadata_df$plate_key %in% unique_samples]
  return(length(unique(matched_participants)))
}

# Filter for canonical variants -------------------------------------------------
# Stats of variant filtering
stats_variants_filter <- data.frame(metric = character(), variant_count = numeric(), 
                                    participant_count = numeric(), stringsAsFactors = FALSE
)

stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_raw_variants", variant_count = nrow(variant_table),
                                participant_count = count_unique_participants(variant_table, participant_metadata)
                              ))

# Filter for canonical variants
variant_filtered_table <- variant_table %>%
  filter(CANONICAL_annotation == "YES")

# Add canonical variants count to stats
stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_canonical_variants", 
                                      variant_count = nrow(variant_filtered_table),
                                      participant_count = count_unique_participants(variant_filtered_table, participant_metadata)
                              ))


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
  category_variants <- variant_filtered_table[variant_filtered_table$clinvar_category == category, ]
  stats_variants_filter <- rbind(stats_variants_filter, 
                                data.frame(metric = paste0("total_canonical_", category), 
                                           variant_count = as.numeric(clinvar_counts[category]),
                                           participant_count = count_unique_participants(category_variants, participant_metadata)))
}

# Filter out benign variants
variant_filtered_table <- variant_filtered_table[variant_filtered_table$clinvar_category != "benign", ]

# Add count after benign filtering
stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_canonical_nonBenign", 
                                         variant_count = nrow(variant_filtered_table),
                                         participant_count = count_unique_participants(variant_filtered_table, participant_metadata)
                              ))


# Filter by MAX_AF annotation --------------------------------------------------
# Keep variants with MAX_AF = "-" or < 0.001
variant_filtered_table <- variant_filtered_table %>%
  filter(!is.na(MAX_AF_annotation) & (MAX_AF_annotation == "-" | as.numeric(MAX_AF_annotation) < 0.001))


# Add count after MAX_AF filtering
stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_canonical_nonBenign_MAF<0.001", 
                                         variant_count = nrow(variant_filtered_table),
                                         participant_count = count_unique_participants(variant_filtered_table, participant_metadata)
                              ))


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
