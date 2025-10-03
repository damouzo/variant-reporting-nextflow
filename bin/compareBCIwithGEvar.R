#!/usr/bin/env Rscript
# Script: compareBCIwithGEvar

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
bci_file    <- args[1]  # "/data/BCI-KRP/projects/RR771-DHX34/results/DDX41/DDX41_bci_patients_variants.rds"
ge_file     <- args[2]  # "/data/BCI-KRP/projects/RR771-DHX34/results/DDX41/DDX41_small_variants.rds"
gene_name   <- args[3]  # "DDX41"


# Libraries  -------------------------------------------------------------------
library(tidyverse)


# Settings ----------------------------------------------------------------------
set.seed(23)


# Load gene annotated variant file ---------------------------------------------
bci_var_table <- readRDS(bci_file)
ge_variant_table <- readRDS(ge_file)
ge_variant_table <- ge_variant_table %>% filter(CANONICAL_annotation == "YES")

# Function definitions ----------------------------------------------------------
# Function to split and compare existing variation names
compare_existing_variations <- function(bci_names, ge_names) {
  if (is.na(bci_names) || is.na(ge_names)) return(FALSE) # Handle NA values

  # Split by comma and trim whitespace
  bci_split <- trimws(strsplit(bci_names, ",")[[1]])
  ge_split <- trimws(strsplit(ge_names, ",")[[1]])
  
  # Check if any name from bci appears in ge
  any(bci_split %in% ge_split)
}

# Function to extract chromosome and position from Location column
extract_location <- function(location_str) {
  if (is.na(location_str)) return(list(chr = NA, pos = NA))
  
  parts <- strsplit(location_str, ":")[[1]]
  if (length(parts) != 2) return(list(chr = NA, pos = NA))
  
  list(chr = parts[1], pos = as.numeric(parts[2]))
}

# Compare variants by existing names ------------------------------------------------
cat("\nComparing variants by existing variation names...\n")
matches_by_name <- data.frame()

if (!is.null(bci_var_table$Existing_variation) && !is.null(ge_variant_table$Existing_variation_annotation)) {
    # Create lookup for valid BCI variants
  bci_valid <- which(!is.na(bci_var_table$Existing_variation))
  
  if (length(bci_valid) > 0) {
    for (i in bci_valid) {
      bci_existing <- bci_var_table$Existing_variation[i]
      
      # Vectorized comparison across GE table
      ge_matches <- sapply(ge_variant_table$Existing_variation_annotation, 
                          function(x) compare_existing_variations(bci_existing, x))
      ge_match_indices <- which(ge_matches)
      
      if (length(ge_match_indices) > 0) {
        match_rows <- data.frame(
          bci_row = rep(i, length(ge_match_indices)),
          ge_row = ge_match_indices,
          match_type = "existing_variation"
        )
        matches_by_name <- rbind(matches_by_name, match_rows)
      }
    }
  }
}

# Compare variants by genomic location ------------------------------------------------
cat("Comparing variants by genomic location...\n")
matches_by_location <- data.frame()

# Extract all BCI locations at once
bci_locations <- lapply(bci_var_table$Location, extract_location)
bci_chrs <- sapply(bci_locations, function(x) x$chr)
bci_pos <- sapply(bci_locations, function(x) x$pos)

# Create lookup for valid BCI variants
bci_valid <- which(!is.na(bci_chrs) & !is.na(bci_pos))

if (length(bci_valid) > 0) {
  # Extract GE chromosomes and positions
  ge_chrs <- as.character(ge_variant_table$CHROM_variant)
  ge_pos <- ge_variant_table$POS_variant
  
  for (i in bci_valid) {
    # Find matching positions vectorized
    chr_matches <- which(bci_chrs[i] == ge_chrs & bci_pos[i] == ge_pos & 
                        !is.na(ge_chrs) & !is.na(ge_pos))
    
    if (length(chr_matches) > 0) {
      match_rows <- data.frame(
        bci_row = rep(i, length(chr_matches)),
        ge_row = chr_matches,
        match_type = "genomic_location"
      )
      matches_by_location <- rbind(matches_by_location, match_rows)
    }
  }
}

# Combine all matches --------------------------------------------------------------
cat("Combining all matches...\n")
all_matches <- rbind(matches_by_name, matches_by_location)

# Remove duplicates (same variant pair found by both methods)
if (nrow(all_matches) > 0) {
  all_matches <- all_matches[!duplicated(all_matches[c("bci_row", "ge_row")]), ]
}

# Create comparison table
if (nrow(all_matches) > 0) {
  # Select columns of interest for BCI and GE
  bci_selected <- bci_var_table[all_matches$bci_row, c("Variant", "Location", "Allele", "Existing_variation", "Patient_REF", "Patient_Centre_REF","CLIN_SIG")]
  ge_selected <- ge_variant_table[all_matches$ge_row, c("ID_variant", "CHROM_variant", "POS_variant", "REF_variant", "ALT_variant", "Existing_variation_annotation", "IMPACT_annotation", "CLIN_SIG_annotation","Het_samples", "Hom_samples", "Hemi_samples")]
  
  # Combine tables
  compared_var_full <- cbind(bci_selected, ge_selected,
                             match_type = all_matches$match_type)
  
  # Remove duplicate variants
  compared_var <- compared_var_full[!duplicated(compared_var_full[c("Variant", "ID_variant")]), ]
  rownames(compared_var) <- NULL # Reset row names
  
  cat(sprintf("Found %d unique matching variants between BCI and GE tables\n", nrow(compared_var)))
  
} else {
  cat("No matching variants found between BCI and GE tables\n")
  
  # Create empty comparison table with appropriate structure
  compared_var <- data.frame(
    Variant = character(0),
    Location = character(0),
    Allele = character(0),
    Existing_variation = character(0),
    Patient_REF = character(0),
    Patient_Centre_REF = character(0),
    CLIN_SIG = character(0),
    ID_variant = character(0),
    CHROM_variant = character(0),
    POS_variant = numeric(0),
    REF_variant = character(0),
    ALT_variant = character(0),
    Existing_variation_annotation = character(0),
    IMPACT_annotation = character(0),
    CLIN_SIG_annotation = character(0),
    Het_samples = character(0),
    Hom_samples = character(0),
    Hemi_samples = character(0),
    match_type = character(0)
  )
}


# Print summary ----------------------------------------------------------------
cat("\nGlimpse of table data of Variants \n")
glimpse(compared_var)

# Count summary of variant and samples 
total_samples <- 0
if (nrow(compared_var) > 0) {
  all_sample_strings <- c(compared_var$Het_samples, compared_var$Hom_samples, compared_var$Hemi_samples)
  all_sample_strings <- all_sample_strings[!is.na(all_sample_strings) & all_sample_strings != "" & all_sample_strings != "-"]
  all_samples <- unique(unlist(strsplit(paste(all_sample_strings, collapse = ","), ",")))
  total_samples <- length(trimws(all_samples))
}

cat("\n=== VARIANT COMPARISON STATISTICS ===\n")
cat("BCI table unique variants:", length(unique(bci_var_table$Variant)), "\n")
cat("GE table unique variants:", length(unique(ge_variant_table$ID_variant)), "\n") 
cat("Total matched variants:", nrow(compared_var), "\n")
cat("Total samples of matched variants:", total_samples, "\n")

# Save variant table  --------------------------------------------------------------
# TSV
out_file_tsv <- paste0(gene_name, "_bci_patients_compared_ge_small_variants.tsv")
write.table(compared_var, file = out_file_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# RDS
out_file_rds <- paste0(gene_name, "_bci_patients_compared_ge_small_variants.rds")
saveRDS(compared_var, file = out_file_rds)

# Print output messages
cat("\nOutput files written:")
cat("\n- TSV file:", out_file_tsv)
cat("\n- RDS file:", out_file_rds, "\n")
