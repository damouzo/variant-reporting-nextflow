#!/usr/bin/env Rscript
# Script: filterStructVar

# Note: The input structural-variant data have already undergone an initial filtering
#       step as part of the GE Structural Variant workflow. The filtering applied in
#       this script is intentionally minimal.
#       Until additional variant- or sample-level information is available, this
#       minimal filter provides sufficient quality for downstream summary and
#       exploratory analyses.

#args <- c("/mnt/c/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/data/raw_data/WGS_Variants/structVar/GRCh38_DDX41_annotated_variants.tsv", "DDX41", 
#"/mnt/c/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/results/DDX41/DDX41_structural_variants_participantMetadata.rds", 
#"filter_older60", "--age_filter", "older_than 60") # # "filter_basic", "filter_onlyHaem", "filter_rmNeuro", etc.

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: filterStructVar.R <clean_structvar_file> <gene_name> <part_metadata_file> <filter_type> [custom_filter_args...]")
}

clean_structvar_file   <- args[1]
gene_name              <- args[2]
part_metadata_file     <- args[3]
filter_type            <- args[4]

# Parse additional custom filter arguments
custom_filter_args <- list()
if (length(args) > 4) {
  additional_args <- args[5:length(args)]
  # Parse arguments in format --arg_name value
  i <- 1
  while (i < length(additional_args)) {
    if (startsWith(additional_args[i], "--")) {
      arg_name <- sub("^--", "", additional_args[i])
      if (i < length(additional_args)) {
        custom_filter_args[[arg_name]] <- additional_args[i + 1]
        i <- i + 2
      } else {
        custom_filter_args[[arg_name]] <- TRUE
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
}

# Libraries  -------------------------------------------------------------------
library(tidyverse)


# Settings ----------------------------------------------------------------------
set.seed(23)


# Load gene annotated variant file ---------------------------------------------
variant_table <- readRDS(clean_structvar_file)
participant_metadata <- readRDS(part_metadata_file)

# Print filter information
cat("STRUCTURAL VARIANT FILTER INFORMATION:\n")
cat("Gene:", gene_name, "\n")
cat("Filter type:", filter_type, "\n")
if (length(custom_filter_args) > 0) {
  cat("Custom filter arguments:\n")
  for (arg_name in names(custom_filter_args)) {
    cat("  -", arg_name, ":", custom_filter_args[[arg_name]], "\n")
  }
}


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

# Custom filtering functions for structural variants ---------------------------
apply_custom_filters <- function(variant_df, metadata_df, custom_args) {
  cat("Applying custom filters to structural variants...\n")
  
  if (length(custom_args) == 0) {
    cat("No custom filters to apply.\n")
    return(variant_df)
  }
  
  # Get all samples from variants to filter by participant metadata
  all_samples <- c()
  samples_in_col <- unlist(strsplit(variant_df$SAMPLE, ","))
  samples_in_col <- trimws(samples_in_col)
  samples_in_col <- samples_in_col[samples_in_col != "NA" & samples_in_col != ""]
  all_samples <- c(all_samples, samples_in_col)
  
  # Get participant IDs for all samples
  variant_participants <- unique(metadata_df$participant_id[metadata_df$plate_key %in% all_samples])
  
  # Apply disease group filters
  if ("disease_group_filter" %in% names(custom_args)) {
    filter_value <- custom_args[["disease_group_filter"]]
    cat("  Applying disease group filter:", filter_value, "\n")
    
    # Get filtered participants based on disease group
    if (filter_value == "haematological_disease") {
      # Include only participants with haematological diseases
      disease_term <- "Haematological and immunological disorders"
      filtered_participants <- metadata_df$participant_id[grepl(disease_term, metadata_df$normalised_disease_group)]
      
    } else if (filter_value == "exclude_neurology") {
      # Exclude participants with neurological diseases
      disease_term <- "Neurology and neurodevelopmental disorders"
      filtered_participants <- metadata_df$participant_id[!metadata_df$normalised_disease_group %in% disease_term]
      
    } else {
      stop("ERROR: The filter_type is not defined in the script")
    }
    
    # Filter variants to only include those from filtered participants
    variant_df <- filter_struct_variants_by_participants(variant_df, metadata_df, filtered_participants)
    cat("    After disease group filtering:", nrow(variant_df), "variants remain\n")
  }
  
  # Apply age filters
  if ("age_filter" %in% names(custom_args)) {
    filter_value <- custom_args[["age_filter"]]
    cat("Applying age filter:", filter_value, "\n")
    
    if (filter_value == "older_than_60") {
      filtered_participants <- metadata_df$participant_id[
        metadata_df$rare_disease_diagnosis_age >= 60 | metadata_df$cancer_diagnosis_age >= 60
      ]
    } else {
      stop("ERROR: The filter_type is not defined in the script")
    }
    
    # Filter variants
    variant_df <- filter_struct_variants_by_participants(variant_df, metadata_df, filtered_participants)
    cat("    After age filtering:", nrow(variant_df), "variants remain\n")
  }
  
  return(variant_df)
}

# Helper function to filter structural variants by participant list
filter_struct_variants_by_participants <- function(variant_df, metadata_df, kept_participants) {
  # Get plate_keys for kept participants
  kept_plate_keys <- metadata_df$plate_key[metadata_df$participant_id %in% kept_participants]
  
  # For each variant, check if it has at least one sample from kept participants
  keep_row <- sapply(1:nrow(variant_df), function(i) {
    # Process SAMPLE column for structural variants
    if ("SAMPLE" %in% colnames(variant_df) && !is.na(variant_df[i, "SAMPLE"]) && variant_df[i, "SAMPLE"] != "NA") {
      samples <- unlist(strsplit(as.character(variant_df[i, "SAMPLE"]), ","))
      samples <- trimws(samples)
      samples <- samples[samples != "NA" & samples != "" & !is.na(samples)]
      
      # Return TRUE if at least one sample is from kept participants
      return(any(samples %in% kept_plate_keys))
    }
    return(FALSE)
  })
  
  # Return only rows where keep_row is TRUE
  return(variant_df[keep_row, ])
}


# Apply STANDARD filtering FIRST ----------------------------------------------
cat("\n=== APPLYING STANDARD FILTERING FOR STRUCTURAL VARIANTS (ALL FILTER TYPES) ===\n")

# Stats of variant filtering
stats_variants_filter <- data.frame(metric = character(), variant_count = numeric(), 
                                    participant_count = numeric(), stringsAsFactors = FALSE)

stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_raw_variants", variant_count = nrow(variant_table),
                                participant_count = count_unique_participants(variant_table, participant_metadata)
                              ))

# Filter for PASS variants ---------
variant_filtered_table <- variant_table[variant_table$FILTER == "PASS", ]

stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_filter_PASS_variants", 
                                      variant_count = nrow(variant_filtered_table),
                                      participant_count = count_unique_participants(variant_filtered_table, participant_metadata)))


# Filter for QUAL_SVonly >= 250 (keep NAs) -----------------
# Reference: https://www.nature.com/articles/s41467-020-16481-5
if ("QUAL_SVonly" %in% colnames(variant_filtered_table)) {
  variant_filtered_table <- variant_filtered_table[is.na(variant_filtered_table$QUAL_SVonly) | 
                                                   variant_filtered_table$QUAL_SVonly >= 250, ]
  
  stats_variants_filter <- rbind(stats_variants_filter, 
                                data.frame(metric = "total_PASS_QUAL250_variants", 
                                           variant_count = nrow(variant_filtered_table),
                                           participant_count = count_unique_participants(variant_filtered_table, participant_metadata)))
}

# Filter for majority gene symbol -----------------
if ("SYMBOL_annotation" %in% colnames(variant_filtered_table)) {
  # Count occurrences of each symbol to keep only majority symbol
  symbol_counts <- table(variant_filtered_table$SYMBOL_annotation)
  majority_symbol <- names(symbol_counts)[which.max(symbol_counts)]
  
  variant_filtered_table <- variant_filtered_table[variant_filtered_table$SYMBOL_annotation == majority_symbol, ]
  
  stats_variants_filter <- rbind(stats_variants_filter, 
                                data.frame(metric = paste0("total_PASS_", majority_symbol, "_named"), 
                                           variant_count = nrow(variant_filtered_table),
                                           participant_count = count_unique_participants(variant_filtered_table, participant_metadata)))
}


# Detailed stats by variant origin and SVTYPE ----------------------------------
# Separate germline and somatic variants
germline_variants <- variant_filtered_table[variant_filtered_table$variant_origin == "germline", ]
somatic_variants <- variant_filtered_table[variant_filtered_table$variant_origin == "somatic", ]

# Add overall germline/somatic counts - same for all filter types
stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_germline_variants", 
                                         variant_count = nrow(germline_variants),
                                         participant_count = count_unique_participants(germline_variants, participant_metadata)))

stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_somatic_variants", 
                                         variant_count = nrow(somatic_variants),
                                         participant_count = count_unique_participants(somatic_variants, participant_metadata)))

cat("High-quality structural variants:", nrow(germline_variants), "germline,", nrow(somatic_variants), "somatic\n")

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

# Store high-quality variants for custom filtering
custom_filtering_variants <- variant_filtered_table

# Apply CUSTOM filtering AFTER standard filtering (if not basic filter) -------
if (filter_type != "filter_basic") {
  cat("\n=== APPLYING CUSTOM FILTERING FOR STRUCTURAL VARIANTS", filter_type, "===\n")
  
  # Apply custom filters to the high-quality variants
  variant_filtered_table <- apply_custom_filters(custom_filtering_variants, participant_metadata, custom_filter_args)
  
  # Add custom filter metric with renamed format
  filter_name <- sub("^filter_", "", filter_type)  # Remove "filter_" prefix
  custom_metric <- paste0("total_", filter_name, "_variants")
  
  stats_variants_filter <- rbind(stats_variants_filter, 
                                data.frame(metric = custom_metric, 
                                           variant_count = nrow(variant_filtered_table),
                                           participant_count = count_unique_participants(variant_filtered_table, participant_metadata)
                                ))
  
  # Add germline and somatic breakdown after custom filtering
  germline_post_custom <- variant_filtered_table[variant_filtered_table$variant_origin == "germline", ]
  somatic_post_custom <- variant_filtered_table[variant_filtered_table$variant_origin == "somatic", ]
  
  stats_variants_filter <- rbind(
    stats_variants_filter,
    data.frame(
      metric = paste0(filter_name, "_germline"),
      variant_count = nrow(germline_post_custom),
      participant_count = count_unique_participants(germline_post_custom, participant_metadata)
    ),
    data.frame(
      metric = paste0(filter_name, "_somatic"),
      variant_count = nrow(somatic_post_custom),
      participant_count = count_unique_participants(somatic_post_custom, participant_metadata)
    )
  )
  
  # Add breakdown by SVTYPE for germline variants after custom filtering
  if (nrow(germline_post_custom) > 0) {
    germline_post_svtype_counts <- table(germline_post_custom$INFO_SVTYPE)
    for (svtype in names(germline_post_svtype_counts)) {
      svtype_variants <- germline_post_custom[germline_post_custom$INFO_SVTYPE == svtype, ]
      stats_variants_filter <- rbind(stats_variants_filter, 
                                    data.frame(metric = paste0(filter_name, "_germline_", svtype), 
                                               variant_count = as.numeric(germline_post_svtype_counts[svtype]),
                                               participant_count = count_unique_participants(svtype_variants, participant_metadata)))
    }
  }
  
  # Add breakdown by SVTYPE for somatic variants after custom filtering
  if (nrow(somatic_post_custom) > 0) {
    somatic_post_svtype_counts <- table(somatic_post_custom$INFO_SVTYPE)
    for (svtype in names(somatic_post_svtype_counts)) {
      svtype_variants <- somatic_post_custom[somatic_post_custom$INFO_SVTYPE == svtype, ]
      stats_variants_filter <- rbind(stats_variants_filter, 
                                    data.frame(metric = paste0(filter_name, "_somatic_", svtype), 
                                               variant_count = as.numeric(somatic_post_svtype_counts[svtype]),
                                               participant_count = count_unique_participants(svtype_variants, participant_metadata)))
    }
  }
} else {
  cat("\n=== BASIC FILTERING - No custom filters applied ===\n")
  # For basic filter, the final result is the high-quality variants
  variant_filtered_table <- custom_filtering_variants
}


# Print summary ----------------------------------------------------------------
cat("\nGlimpse of table data of Structural Variants \n")
glimpse(variant_filtered_table)


# Save variant table  --------------------------------------------------------------
# Save stats file with filter type in filename
stats_file_csv <- paste0(gene_name, "_structural_variants_", filter_type, "_stats.csv")
write.csv(stats_variants_filter, file = stats_file_csv, row.names = FALSE)

# TSV with filter type in filename
out_file_tsv <- paste0(gene_name, "_structural_variants_", filter_type, "_filtered.tsv")
write.table(variant_filtered_table, file = out_file_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# RDS with filter type in filename
out_file_rds <- paste0(gene_name, "_structural_variants_", filter_type, "_filtered.rds")
saveRDS(variant_filtered_table, file = out_file_rds)

# Print output messages
cat("STRUCTURAL VARIANT FILTERING COMPLETED - ", filter_type, "\n")
cat("Gene:", gene_name, "\n")
cat("Final variant count:", nrow(variant_filtered_table), "\n")
cat("Final participant count:", count_unique_participants(variant_filtered_table, participant_metadata), "\n")
cat("\nOutput files written:")
cat("\n- Stats file:", stats_file_csv)
cat("\n- TSV file:", out_file_tsv)
cat("\n- RDS file:", out_file_rds, "\n")
