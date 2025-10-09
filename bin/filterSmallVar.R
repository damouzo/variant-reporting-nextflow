#!/usr/bin/env Rscript
# Script: filterSmallVar

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: filterSmallVar.R <clean_smallvar_file> <gene_name> <part_metadata_file> <filter_type> [custom_filter_args...]")
}

clean_smallvar_file   <- args[1]  # "/mnt/c/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/results/DHX34/DHX34_small_variants.rds
gene_name             <- args[2]  # "DHX34"
part_metadata_file    <- args[3]  # "/mnt/c/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/results/DHX34/DHX34_structural_variants_participantMetadata.rds"
filter_type           <- args[4]  # "filter_basic", "filter_onlyHaem", "filter_rmNeuro", etc.

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
variant_table <- readRDS(clean_smallvar_file)
participant_metadata <- readRDS(part_metadata_file)

# Print filter information
cat("\n" %>% str_dup(50), "\n")
cat("FILTER INFORMATION:\n")
cat("Gene:", gene_name, "\n")
cat("Filter type:", filter_type, "\n")
if (length(custom_filter_args) > 0) {
  cat("Custom filter arguments:\n")
  for (arg_name in names(custom_filter_args)) {
    cat("  -", arg_name, ":", custom_filter_args[[arg_name]], "\n")
  }
}
cat("" %>% str_dup(50), "\n")


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

# Custom filtering functions ---------------------------------------------------
apply_custom_filters <- function(variant_df, metadata_df, custom_args) {
  cat("Applying custom filters...\n")
  
  if (length(custom_args) == 0) {
    cat("No custom filters to apply.\n")
    return(variant_df)
  }
  
  # Get all samples from variants to filter by participant metadata
  all_samples <- c()
  for (col in c("Het_samples", "Hom_samples", "Hemi_samples")) {
    if (col %in% colnames(variant_df)) {
      samples_in_col <- unlist(strsplit(variant_df[[col]], ","))
      samples_in_col <- trimws(samples_in_col)
      samples_in_col <- samples_in_col[samples_in_col != "NA" & samples_in_col != ""]
      all_samples <- c(all_samples, samples_in_col)
    }
  }
  
  # Get participant IDs for all samples
  variant_participants <- unique(metadata_df$participant_id[metadata_df$plate_key %in% all_samples])
  
  # Apply disease group filters
  if ("disease_group_filter" %in% names(custom_args)) {
    filter_value <- custom_args[["disease_group_filter"]]
    cat("  Applying disease group filter:", filter_value, "\n")
    
    # Get filtered participants based on disease group
    if (filter_value == "haematological_disease" || filter_value == "haematological_malignancies") {
      # Keep only participants with hematological diseases
      haem_diseases <- c("Acute leukaemia", "Chronic leukaemia", "Myelodysplastic syndrome", 
                        "Myeloproliferative neoplasm", "Lymphoma", "Multiple myeloma", 
                        "Aplastic anaemia", "Primary immunodeficiency", "Haemoglobinopathy",
                        "Primary haemostasis and coagulation disorders", "Other haematological disorders")
      
      if ("disease_group" %in% colnames(metadata_df)) {
        filtered_participants <- metadata_df$participant_id[metadata_df$disease_group %in% haem_diseases]
      } else {
        cat("  Warning: 'disease_group' column not found in metadata. Skipping filter.\n")
        filtered_participants <- variant_participants
      }
      
    } else if (filter_value == "exclude_neurology") {
      # Exclude participants with neurological diseases
      neuro_diseases <- c("Neurology", "Neuromuscular disorders", "Inherited neurological disorders", 
                         "Epilepsy", "Intellectual disability", "Autism spectrum disorder",
                         "Developmental disorders")
      
      if ("disease_group" %in% colnames(metadata_df)) {
        filtered_participants <- metadata_df$participant_id[!metadata_df$disease_group %in% neuro_diseases]
      } else {
        cat("  Warning: 'disease_group' column not found in metadata. Skipping filter.\n")
        filtered_participants <- variant_participants
      }
    } else {
      cat("  Warning: Unknown disease group filter:", filter_value, ". Skipping.\n")
      filtered_participants <- variant_participants
    }
    
    # Filter variants to only include those from filtered participants
    variant_df <- filter_variants_by_participants(variant_df, metadata_df, filtered_participants)
    cat("    After disease group filtering:", nrow(variant_df), "variants remain\n")
  }
  
  # Apply age filters
  if ("age_filter" %in% names(custom_args)) {
    filter_value <- custom_args[["age_filter"]]
    cat("  Applying age filter:", filter_value, "\n")
    
    if (filter_value == "older_than_60") {
      if ("age" %in% colnames(metadata_df)) {
        filtered_participants <- metadata_df$participant_id[!is.na(metadata_df$age) & metadata_df$age > 60]
      } else {
        cat("  Warning: 'age' column not found in metadata. Skipping filter.\n")
        filtered_participants <- variant_participants
      }
      
      # Filter variants
      variant_df <- filter_variants_by_participants(variant_df, metadata_df, filtered_participants)
      cat("    After age filtering:", nrow(variant_df), "variants remain\n")
    }
  }
  
  return(variant_df)
}

# Helper function to filter variants by participant list
filter_variants_by_participants <- function(variant_df, metadata_df, kept_participants) {
  # Get plate_keys for kept participants
  kept_plate_keys <- metadata_df$plate_key[metadata_df$participant_id %in% kept_participants]
  
  # Filter each row to only include samples from kept participants
  variant_df_filtered <- variant_df[0, ]  # Empty dataframe with same structure
  
  for (i in 1:nrow(variant_df)) {
    row <- variant_df[i, ]
    
    # Process each sample type column
    for (col in c("Het_samples", "Hom_samples", "Hemi_samples")) {
      if (col %in% colnames(row) && !is.na(row[[col]]) && row[[col]] != "NA") {
        samples <- unlist(strsplit(as.character(row[[col]]), ","))
        samples <- trimws(samples)
        samples <- samples[samples != "NA" & samples != ""]
        
        # Keep only samples that are in kept_plate_keys
        filtered_samples <- samples[samples %in% kept_plate_keys]
        
        # Update the column
        if (length(filtered_samples) > 0) {
          row[[col]] <- paste(filtered_samples, collapse = ",")
        } else {
          row[[col]] <- "NA"
        }
      }
    }
    
    # Only keep the row if at least one sample column has valid samples
    has_valid_samples <- FALSE
    for (col in c("Het_samples", "Hom_samples", "Hemi_samples")) {
      if (col %in% colnames(row) && !is.na(row[[col]]) && row[[col]] != "NA") {
        has_valid_samples <- TRUE
        break
      }
    }
    
    if (has_valid_samples) {
      variant_df_filtered <- rbind(variant_df_filtered, row)
    }
  }
  
  return(variant_df_filtered)
}

# Apply STANDARD filtering FIRST (for all filter types) -----------------------
# This ensures all filter types start from the same high-quality baseline
cat("\n=== APPLYING STANDARD FILTERING (ALL FILTER TYPES) ===\n")

# Stats of variant filtering - all filter types use same baseline for standard steps
stats_variants_filter <- data.frame(metric = character(), variant_count = numeric(), 
                                    participant_count = numeric(), stringsAsFactors = FALSE
)

# Initial raw count - same for all filter types for comparability
stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_raw_variants", variant_count = nrow(variant_table),
                                participant_count = count_unique_participants(variant_table, participant_metadata)
                              ))

cat("Starting with", nrow(variant_table), "raw variants\n")

# Filter for canonical variants
variant_filtered_table <- variant_table %>%
  filter(CANONICAL_annotation == "YES")

cat("After canonical filtering:", nrow(variant_filtered_table), "variants remain\n")

# Add canonical variants count to stats - same for all filter types
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

# Count each category and add to stats - same for all filter types
clinvar_counts <- table(variant_filtered_table$clinvar_category)
for (category in names(clinvar_counts)) {
  category_variants <- variant_filtered_table[variant_filtered_table$clinvar_category == category, ]
  category_metric <- paste0("total_canonical_", category)
  stats_variants_filter <- rbind(stats_variants_filter, 
                                data.frame(metric = category_metric, 
                                           variant_count = as.numeric(clinvar_counts[category]),
                                           participant_count = count_unique_participants(category_variants, participant_metadata)))
}

# Filter out benign variants
variant_filtered_table <- variant_filtered_table[variant_filtered_table$clinvar_category != "benign", ]

cat("After removing benign variants:", nrow(variant_filtered_table), "variants remain\n")

# Add count after benign filtering - same for all filter types
stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_canonical_nonBenign", 
                                         variant_count = nrow(variant_filtered_table),
                                         participant_count = count_unique_participants(variant_filtered_table, participant_metadata)
                              ))


# Filter by MAX_AF annotation --------------------------------------------------
# Keep variants with MAX_AF = "-" or < 0.001
variant_filtered_table <- variant_filtered_table %>%
  filter(!is.na(MAX_AF_annotation) & (MAX_AF_annotation == "-" | as.numeric(MAX_AF_annotation) < 0.001))

cat("After MAF filtering:", nrow(variant_filtered_table), "high-quality variants remain\n")

# Add count after MAX_AF filtering - same for all filter types (this is the high-quality baseline)
stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_canonical_nonBenign_MAF<0.001", 
                                         variant_count = nrow(variant_filtered_table),
                                         participant_count = count_unique_participants(variant_filtered_table, participant_metadata)
                              ))

# Store high-quality variants for custom filtering
high_quality_variants <- variant_filtered_table

# Apply CUSTOM filtering AFTER standard filtering (if not basic filter) -------
if (filter_type != "filter_basic") {
  cat("\n=== APPLYING CUSTOM FILTERING FOR", filter_type, "===\n")
  
  # Apply custom filters to the high-quality variants
  variant_filtered_table <- apply_custom_filters(high_quality_variants, participant_metadata, custom_filter_args)
  
  cat("After custom filtering:", nrow(variant_filtered_table), "variants remain\n")
  
  # Add custom filter metric with specific prefix
  custom_metric <- paste0(filter_type, "_after_custom_filtering")
  stats_variants_filter <- rbind(stats_variants_filter, 
                                data.frame(metric = custom_metric, 
                                           variant_count = nrow(variant_filtered_table),
                                           participant_count = count_unique_participants(variant_filtered_table, participant_metadata)
                                ))
} else {
  cat("\n=== BASIC FILTERING - No custom filters applied ===\n")
  # For basic filter, the final result is the high-quality variants
  variant_filtered_table <- high_quality_variants
}


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
cat("\n" %>% str_dup(50), "\n")
cat("FILTERING COMPLETED - ", filter_type, "\n")
cat("Gene:", gene_name, "\n")
cat("Final variant count:", nrow(variant_filtered_table), "\n")
cat("Final participant count:", count_unique_participants(variant_filtered_table, participant_metadata), "\n")
cat("\nOutput files written:")
cat("\n- Stats file:", stats_file_csv)
cat("\n- TSV file:", out_file_tsv)
cat("\n- RDS file:", out_file_rds, "\n")
cat("" %>% str_dup(50), "\n")
