#!/usr/bin/env Rscript
# Script: filterStructVar

# Note: The input structural-variant data have already undergone an initial filtering
#       step as part of the GE Structural Variant workflow. The filtering applied in
#       this script is intentionally minimal.
#       Until additional variant- or sample-level information is available, this
#       minimal filter provides sufficient quality for downstream summary and
#       exploratory analyses.

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
clean_structvar_file   <- args[1]  # "C:/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/data/raw_data/WGS_Variants/structVar/GRCh38_DDX41_annotated_variants.tsv"
gene_name              <- args[2]        # "DDX41"
part_metadata_file     <- args[3]   # TSV file with participant metadata

# Libraries  -------------------------------------------------------------------
library(tidyverse)


# Settings ----------------------------------------------------------------------
set.seed(23)


# Load gene annotated variant file ---------------------------------------------
variant_table <- readRDS(clean_structvar_file)
participant_metadata <- read.table(part_metadata_file, sep = "\t", header = TRUE, 
                                  stringsAsFactors = FALSE, quote = "")


# Function to count unique participants ----------------------------------------
count_unique_participants <- function(variant_df, metadata_df) {
  # Extract all samples from SAMPLE column
  all_samples <- c()
  
  if ("SAMPLE" %in% colnames(variant_df)) {
    samples_in_col <- unlist(strsplit(variant_df$SAMPLE, ","))
    samples_in_col <- trimws(samples_in_col)
    samples_in_col <- samples_in_col[samples_in_col != "NA" & samples_in_col != ""]
    all_samples <- c(all_samples, samples_in_col)
  }
  
  # Get unique participant_ids
  unique_samples <- unique(all_samples)
  matched_participants <- metadata_df$participant_id[metadata_df$plate_key %in% unique_samples]
  return(length(unique(matched_participants)))
}


# Filter for PASS variants and majority gene symbol ----------------------------
# Stats of variant filtering
stats_variants_filter <- data.frame(metric = character(), variant_count = numeric(), 
                                    participant_count = numeric(), stringsAsFactors = FALSE)

stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_raw_variants", variant_count = nrow(variant_table),
                                participant_count = count_unique_participants(variant_table, participant_metadata)))

# Filter for PASS variants
variant_filtered_table <- variant_table[variant_table$FILTER == "PASS", ]

stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_filter_PASS_variants", 
                                      variant_count = nrow(variant_filtered_table),
                                      participant_count = count_unique_participants(variant_filtered_table, participant_metadata)))

# Filter for QUAL_SVonly >= 250 (keep NAs)
# Reference: https://www.nature.com/articles/s41467-020-16481-5
if ("QUAL_SVonly" %in% colnames(variant_filtered_table)) {
  variant_filtered_table <- variant_filtered_table[is.na(variant_filtered_table$QUAL_SVonly) | 
                                                   variant_filtered_table$QUAL_SVonly >= 250, ]
  
  stats_variants_filter <- rbind(stats_variants_filter, 
                                data.frame(metric = "total_PASS_QUAL250_variants", 
                                           variant_count = nrow(variant_filtered_table),
                                           participant_count = count_unique_participants(variant_filtered_table, participant_metadata)))
}

# Filter for majority gene symbol
if ("SYMBOL_annotation" %in% colnames(variant_filtered_table)) {
  # Count occurrences of each symbol
  symbol_counts <- table(variant_filtered_table$SYMBOL_annotation)
  majority_symbol <- names(symbol_counts)[which.max(symbol_counts)]
  
  # Filter to keep only majority symbol
  variant_filtered_table <- variant_filtered_table[variant_filtered_table$SYMBOL_annotation == majority_symbol, ]
  
  stats_variants_filter <- rbind(stats_variants_filter, 
                                data.frame(metric = paste0("total_PASS_", majority_symbol), 
                                           variant_count = nrow(variant_filtered_table),
                                           participant_count = count_unique_participants(variant_filtered_table, participant_metadata)))
}


# Detailed stats by variant origin and SVTYPE ----------------------------------
# Separate germline and somatic variants
germline_variants <- variant_filtered_table[variant_filtered_table$variant_origin == "germline", ]
somatic_variants <- variant_filtered_table[variant_filtered_table$variant_origin == "somatic", ]

# Add overall germline/somatic counts
stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_germline_variants", 
                                         variant_count = nrow(germline_variants),
                                         participant_count = count_unique_participants(germline_variants, participant_metadata)))

stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_somatic_variants", 
                                         variant_count = nrow(somatic_variants),
                                         participant_count = count_unique_participants(somatic_variants, participant_metadata)))

# Count by SVTYPE for germline variants
if (nrow(germline_variants) > 0) {
  germline_svtype_counts <- table(germline_variants$INFO_SVTYPE)
  for (svtype in names(germline_svtype_counts)) {
    svtype_variants <- germline_variants[germline_variants$INFO_SVTYPE == svtype, ]
    stats_variants_filter <- rbind(stats_variants_filter, 
                                  data.frame(metric = paste0("germline_", svtype), 
                                             variant_count = as.numeric(germline_svtype_counts[svtype]),
                                             participant_count = count_unique_participants(svtype_variants, participant_metadata)))
  }
}

# Count by SVTYPE for somatic variants
if (nrow(somatic_variants) > 0) {
  somatic_svtype_counts <- table(somatic_variants$INFO_SVTYPE)
  for (svtype in names(somatic_svtype_counts)) {
    svtype_variants <- somatic_variants[somatic_variants$INFO_SVTYPE == svtype, ]
    stats_variants_filter <- rbind(stats_variants_filter, 
                                  data.frame(metric = paste0("somatic_", svtype), 
                                             variant_count = as.numeric(somatic_svtype_counts[svtype]),
                                             participant_count = count_unique_participants(svtype_variants, participant_metadata)))
  }
}


# Print summary ----------------------------------------------------------------
cat("\nGlimpse of table data of Structural Variants \n")
glimpse(variant_filtered_table)


# Save variant table  --------------------------------------------------------------
# Save stats file
stats_file_csv <- paste0(gene_name, "_struct_variants_filtered_stats.csv")
write.csv(stats_variants_filter, file = stats_file_csv, row.names = FALSE)

# TSV
out_file_tsv <- paste0(gene_name, "_struct_variants_filtered.tsv")
write.table(variant_filtered_table, file = out_file_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# RDS
out_file_rds <- paste0(gene_name, "_struct_variants_filtered.rds")
saveRDS(variant_filtered_table, file = out_file_rds)

# Print output messages
cat("\nOutput files written:")
cat("\n- Stats file:", stats_file_csv)
cat("\n- TSV file:", out_file_tsv)
cat("\n- RDS file:", out_file_rds, "\n")
