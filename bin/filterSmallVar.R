#!/usr/bin/env Rscript
# Script: filterSmallVar

#args <- c("/mnt/c/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/results/DHX34/DHX34_small_variants.rds", "DHX34", 
#"/mnt/c/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/results/DHX34/DHX34_structural_variants_participantMetadata.rds", 
#"filter_older60", "--age_filter", "older_than 60")

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: filterSmallVar.R <clean_smallvar_file> <gene_name> <part_metadata_file> <filter_type> [custom_filter_args...]")
}

clean_smallvar_file   <- args[1]
gene_name             <- args[2]
part_metadata_file    <- args[3]
filter_type           <- args[4]

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
cat("FILTER INFORMATION:\n")
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
  # Extract all samples from Het_samples, Hom_samples, and Hemi_samples
  all_samples <- c()
  
  for (col in c("Het_samples", "Hom_samples", "Hemi_samples")) {
    samples_in_col <- unlist(strsplit(variant_df[[col]], ","))
    samples_in_col <- trimws(samples_in_col)
    samples_in_col <- samples_in_col[samples_in_col != "NA" & samples_in_col != ""]
    all_samples <- c(all_samples, samples_in_col)
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
    samples_in_col <- unlist(strsplit(variant_df[[col]], ","))
    samples_in_col <- trimws(samples_in_col)
    samples_in_col <- samples_in_col[samples_in_col != "NA" & samples_in_col != ""]
    all_samples <- c(all_samples, samples_in_col)
  }
  
  # Get participant IDs for all samples
  variant_participants <- unique(metadata_df$participant_id[metadata_df$plate_key %in% all_samples])
  
  # Apply disease group filters
  if ("disease_group_filter" %in% names(custom_args)) {
    filter_value <- custom_args[["disease_group_filter"]]
    cat("Applying disease group filter:", filter_value, "\n")
    
    # Get filtered participants based on disease group
    if (filter_value == "haematological_disease") {
      # Include only participants with haematological diseases
      disease_term <- "Haematological and immunological disorders"
      filtered_participants <- metadata_df$participant_id[grepl(disease_term, metadata_df$normalised_disease_group)]
      
    } else if (filter_value == "exclude_neurology") {
      # Exclude participants with neurological diseases
      disease_term <- "Neurology and neurodevelopmental disorders"
      filtered_participants <- metadata_df$participant_id[!grepl(disease_term, metadata_df$normalised_disease_group)]
      
    } else {
      stop("ERROR: The filter_type is not defined in the script")
    }
    
    # Filter variants to only include those from filtered participants
    variant_df <- filter_variants_by_participants(variant_df, metadata_df, filtered_participants)
    cat("After disease group filtering:", nrow(variant_df), "variants remain\n")
  }
  
  # Apply age filters
  if ("age_filter" %in% names(custom_args)) {
    filter_value <- custom_args[["age_filter"]]
    cat("Applying age filter:", filter_value, "\n")
    
    if (filter_value == "older_than_60" || filter_value == "older_than 60") {
      filtered_participants <- metadata_df$participant_id[
        metadata_df$rare_disease_diagnosis_age >= 60 | metadata_df$cancer_diagnosis_age >= 60
      ]
    } else {
      stop("ERROR: The filter_type is not defined in the script")
    }
    
    # Filter variants
    variant_df <- filter_variants_by_participants(variant_df, metadata_df, filtered_participants)
    cat("After age filtering:", nrow(variant_df), "variants remain\n")
  }
  return(variant_df)
}

# Helper function to filter variants by participant list
filter_variants_by_participants <- function(variant_df, metadata_df, kept_participants) {
  # Get plate_keys for kept participants
  kept_plate_keys <- metadata_df$plate_key[metadata_df$participant_id %in% kept_participants]
  
  # For each variant, check if it has at least one sample from kept participants
  keep_row <- sapply(1:nrow(variant_df), function(i) {
    # Get all samples from this variant
    all_samples <- c()
    for (col in c("Het_samples", "Hom_samples", "Hemi_samples")) {
      if (col %in% colnames(variant_df)) {
        samples <- unlist(strsplit(as.character(variant_df[i, col]), ","))
        samples <- trimws(samples)
        samples <- samples[samples != "NA" & samples != "" & !is.na(samples)]
        all_samples <- c(all_samples, samples)
      }
    }
    # Return TRUE if at least one sample is from kept participants
    return(any(all_samples %in% kept_plate_keys))
  })
  # Return only rows where keep_row is TRUE
  return(variant_df[keep_row, ])
}

# Apply STANDARD filtering FIRST ----------------------------------------------
cat("\n=== APPLYING STANDARD FILTERING (ALL FILTER TYPES) ===\n")

# Stats of variant filtering
stats_variants_filter <- data.frame(metric = character(), variant_count = numeric(), 
                                    participant_count = numeric(), stringsAsFactors = FALSE
)

stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_raw_variants", variant_count = nrow(variant_table),
                                participant_count = count_unique_participants(variant_table, participant_metadata)
                              ))

## Filter for canonical variants -----------
variant_filtered_table <- variant_table %>% filter(CANONICAL_annotation == "YES")

stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_canonical_variants", 
                                      variant_count = nrow(variant_filtered_table),
                                      participant_count = count_unique_participants(variant_filtered_table, participant_metadata)
                              ))


## Filter ClinVar annotations -------
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
  category_metric <- paste0("total_canonical_", category)
  stats_variants_filter <- rbind(stats_variants_filter, 
                                data.frame(metric = category_metric, 
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

# Add count after MAX_AF filtering  (this is the high-quality baseline)
stats_variants_filter <- rbind(stats_variants_filter, 
                              data.frame(metric = "total_canonical_nonBenign_MAF<0.001", 
                                         variant_count = nrow(variant_filtered_table),
                                         participant_count = count_unique_participants(variant_filtered_table, participant_metadata)
                              ))

# Store high-quality variants for custom filtering
custom_filtering_variants <- variant_filtered_table

# Apply CUSTOM filtering AFTER standard filtering (if not basic filter) -------
# Initialize filtered_metadata to be the same as original for basic filter
filtered_metadata <- participant_metadata

if (filter_type != "filter_basic") {
  cat("\n=== APPLYING CUSTOM FILTERING FOR", filter_type, "===\n")
  
  # Apply custom filters to the high-quality variants
  variant_filtered_table <- apply_custom_filters(custom_filtering_variants, participant_metadata, custom_filter_args)
  
  # Also filter the metadata to match the filtering logic
  if ("disease_group_filter" %in% names(custom_filter_args)) {
    filter_value <- custom_filter_args[["disease_group_filter"]]
    
    if (filter_value == "haematological_disease") {
      # Include only participants with haematological diseases
      disease_term <- "Haematological and immunological disorders"
      filtered_metadata <- participant_metadata %>%
        filter(grepl(disease_term, normalised_disease_group))
      
    } else if (filter_value == "exclude_neurology") {
      # Exclude participants with neurological diseases
      disease_term <- "Neurology and neurodevelopmental disorders"
      filtered_metadata <- participant_metadata %>%
        filter(!grepl(disease_term, normalised_disease_group))
    }
    
    cat("Metadata filtered - participants remaining:", nrow(filtered_metadata), "\n")
  }
  
  # Apply age filters to metadata if specified
  if ("age_filter" %in% names(custom_filter_args)) {
    filter_value <- custom_filter_args[["age_filter"]]
    
    if (filter_value == "older_than_60" || filter_value == "older_than 60") {
      filtered_metadata <- filtered_metadata %>%
        filter(rare_disease_diagnosis_age >= 60 | cancer_diagnosis_age >= 60)
      cat("Age filter applied to metadata - participants remaining:", nrow(filtered_metadata), "\n")
    }
  }
  
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
  variant_filtered_table <- custom_filtering_variants
  # For basic filter, metadata remains unchanged
  filtered_metadata <- participant_metadata
}


# Print summary ----------------------------------------------------------------
cat("\nGlimpse of table data of Variants \n")
glimpse(variant_filtered_table)


# Save variant table  --------------------------------------------------------------
# Save stats file with filter type in filename
stats_file_csv <- paste0(gene_name, "_small_variants_", filter_type, "_stats.csv")
write.csv(stats_variants_filter, file = stats_file_csv, row.names = FALSE)

# TSV with filter type in filename
out_file_tsv <- paste0(gene_name, "_small_variants_", filter_type, "_filtered.tsv")
write.table(variant_filtered_table, file = out_file_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# RDS with filter type in filename
out_file_rds <- paste0(gene_name, "_small_variants_", filter_type, "_filtered.rds")
saveRDS(variant_filtered_table, file = out_file_rds)

# IMPORTANT: Save filtered metadata
metadata_output_file <- paste0(gene_name, "_small_variants_", filter_type, "_filtered_metadata.rds")
saveRDS(filtered_metadata, metadata_output_file)

# Print output messages
cat("FILTERING COMPLETED - ", filter_type, "\n")
cat("Gene:", gene_name, "\n")
cat("Final variant count:", nrow(variant_filtered_table), "\n")
cat("Final participant count:", count_unique_participants(variant_filtered_table, participant_metadata), "\n")
cat("\nOutput files written:")
cat("\n- Stats file:", stats_file_csv)
cat("\n- TSV file:", out_file_tsv)
cat("\n- RDS file:", out_file_rds)
cat("\n- Filtered metadata file:", metadata_output_file, "\n")
