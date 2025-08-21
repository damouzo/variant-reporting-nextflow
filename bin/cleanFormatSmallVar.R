#!/usr/bin/env Rscript
# Script: cleanFormatSmallVar

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
annotation_file <- args[1]  # "C:/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/data/raw_data/WGS_Variants/smallVar/GRCh38_DDX41_annotated_variants.tsv"
gene_name <- args[2]        # "DDX41"

# Libraries  -------------------------------------------------------------------
library(tidyverse)


# Settings ----------------------------------------------------------------------
set.seed(23)


# Load gene annotated variant file ---------------------------------------------
variant_table <- read_tsv(annotation_file)


# Clean Data -------------------------------------------------------------------
# Deal with column format issues
columns_to_character <- c("CHROM_variant", "Location_annotation", "Het_samples",
              "Hemi_samples", "Hom_samples")
columns_to_numeric <- c("TSL_annotation")

variant_table <- variant_table %>%
  mutate(across(all_of(columns_to_character), as.character),
     across(all_of(columns_to_numeric), as.numeric))


# Take care of possible ? characters in HGNC_ID
if ("HGNC_ID_annotation" %in% names(variant_table)) {
  if (is.character(variant_table$HGNC_ID_annotation)) {
    variant_table$HGNC_ID_annotation <- sub("^.*?:", "", variant_table$HGNC_ID_annotation)
    variant_table$HGNC_ID_annotation <- as.numeric(variant_table$HGNC_ID_annotation)
  }
}


# Factor with order Impact
variant_table$IMPACT_annotation <- factor(
  variant_table$IMPACT_annotation,
  levels = c("HIGH", "MODERATE", "LOW", "MODIFIER")
)

# Protein Position with Single Number(Taking care of "?-23")
variant_table <- variant_table %>%
  mutate(
    Protein_pos_start = case_when(
      str_detect(Protein_position_annotation, "^\\?-\\d+$") ~ as.numeric(str_extract(Protein_position_annotation, "\\d+$")),
      str_detect(Protein_position_annotation, "^\\d+-\\d+$") ~ as.numeric(str_extract(Protein_position_annotation, "^\\d+")),
      TRUE ~ as.numeric(Protein_position_annotation)
    )
  )

# Add naming Labeling Variants when not available existing var annot
variant_table <- variant_table %>%
  mutate(
    LabelVarPlot = Existing_variation_annotation,
    LabelVarPlot = ifelse(
      LabelVarPlot == "-" | is.na(LabelVarPlot),
      paste0(CHROM_variant, "_", POS_variant, "_", REF_variant, "/", ALT_variant),
      LabelVarPlot
    )
  )


# Print summary ----------------------------------------------------------------
cat("\nGlimpse of table data of Variants \n")
glimpse(variant_table)


# Save variant table  --------------------------------------------------------------
# TSV
out_file_tsv <- paste0(gene_name, "_small_variants.tsv")
write.table(variant_table, file = out_file_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# RDS
out_file_rds <- paste0(gene_name, "_small_variants.rds")
saveRDS(variant_table, file = out_file_rds)

# Print output messages
cat("\nOutput files written:")
cat("\n- TSV file:", out_file_tsv)
cat("\n- RDS file:", out_file_rds, "\n")
