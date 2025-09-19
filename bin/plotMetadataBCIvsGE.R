#!/usr/bin/env Rscript
# Script: plotQCsmallVar

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: plotMetadataBCIvsGE.R <gene_name> <comparison_rds> <bci_clean_rds> <part_metadata_file>")
}

gene_name           <- args[1] # "DDX41"
comparison_file     <- args[2] # "DDX41_compared_bci_ge_small_variants.rds"
bci_clean_file      <- args[3] # "DDX41_bci_patients_variants.rds"
metadata_file       <- args[4] # "DDX41_participant_metadata.tsv"


# Message validation files
cat("Validating input files for gene:", gene_name, "\n")
missing_files <- c()
if (!file.exists(comparison_file)) missing_files <- c(missing_files, paste("Comparison file:", comparison_file))
if (!file.exists(bci_clean_file)) missing_files <- c(missing_files, paste("BCI clean file:", bci_clean_file))
if (!file.exists(metadata_file)) missing_files <- c(missing_files, paste("Metadata file:", metadata_file))

if (length(missing_files) > 0) {
  cat("ERROR: Missing files detected:\n")
  cat(paste(missing_files, collapse = "\n"), "\n")
  stop("Cannot proceed with missing input files")
}

cat("All input files validated successfully\n")

# Libraries  -------------------------------------------------------------------
library(tidyverse)
library(paletteer)
library(patchwork) 


# Settings ----------------------------------------------------------------------
set.seed(23)

# Palettete ---------------------------------------------------------------------
my_pal <- c(paletteer_d("RColorBrewer::Set1"),paletteer_d("RColorBrewer::Set3"))
pastel_colors <- c("#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", "#BAE1FF", "#D5BAFF", 
                   "#FFC8E1", "#C2F0FC","#C5FAD5", "#FFDAC1", "#E2F0CB", "#C3B1E1", "#F8B195", 
                   "#F67280","#F6E2B3", "#B5EAD7","#FFABAB", "#D0F4DE", "#A0CED9", "#FFCBF2")


# Load Data ---------------------------------------------------------------------
# Comparison data
comparison_data <- readRDS(comparison_file)

# BCI clean variants annotated
bci_clean_data <- readRDS(bci_clean_file)

# Metadata from GE
metadata <- read_tsv(metadata_file)



# GE samples with BCI variants --------------------------------------------------
ge_samples_with_bci_var <- comparison_data %>% 
  select(Existing_variation_annotation, Het_samples, Hom_samples, Hemi_samples) %>% 
  mutate(rs_id = str_extract(Existing_variation_annotation, "rs[0-9]+")) %>%
  mutate(across(c(Het_samples, Hom_samples, Hemi_samples), ~na_if(., "NA"))) %>%
  rowwise() %>%
  mutate(all_samples = paste(c(Het_samples, Hom_samples, Hemi_samples)[!is.na(c(Het_samples, Hom_samples, Hemi_samples))], collapse = ",")) %>%
  ungroup() %>%
  separate_rows(all_samples, sep = ",") %>%
  mutate(all_samples = str_trim(all_samples)) %>%
  filter(!is.na(all_samples), all_samples != "") %>% 
  reframe(samples = list(unique(all_samples)), .by = rs_id) %>% 
  { set_names(.$samples, .$rs_id) }
ge_samples_with_bci_var$AllVariantsBCI <- ge_samples_with_bci_var %>% 
  unlist(use.names = FALSE) %>% 
  unique()

# Save in csv files
for (rs_id in names(ge_samples_with_bci_var)) {
  df <- ge_samples_with_bci_var[[rs_id]]
  write_csv(data.frame(df), paste0(gene_name, "_GEsamples_with_BCIvariants_", rs_id, ".csv"), col_names = FALSE)
}



# GE patients metadata with BCI variants ----------------------------------------
ge_patients_with_bci_var <- lapply(ge_samples_with_bci_var, function(samples) {
  metadata %>% filter(plate_key %in% samples)
})

# Save in csv files
for (rs_id in names(ge_patients_with_bci_var)) {
  df <- ge_patients_with_bci_var[[rs_id]]
  write_csv(df, paste0(gene_name, "_GEpatientMetadata_with_BCIvariants_", rs_id, ".csv"))
}


# Plots ------------------------------------------------------------------------
# 1. Distribution of BCI variants in GE patients
vars_of_interest <- c("affection_status", "participant_type", "program", 
                      "participant_karyotype_sex", "normalised_disease_group", 
                      "normalised_specific_disease")

for (rs_id in names(ge_patients_with_bci_var)) {
  df <- ge_patients_with_bci_var[[rs_id]]
  
  plot_list <- lapply(vars_of_interest, function(var) {
    if (!var %in% colnames(df)) return(NULL)
    
    df_count <- df %>%
      filter(!is.na(.data[[var]])) %>%
      count(.data[[var]]) %>%
      mutate(pct = n / sum(n) * 100,
             label = ifelse(pct < 5, "<5%", paste0(round(pct,1), "%")))
    
    # Calculate y-axis limit to accommodate labels
    max_pct <- max(df_count$pct)
    y_limit <- max_pct * 1.15 
    
    ggplot(df_count, aes(x = .data[[var]], y = pct, fill = .data[[var]])) +
      geom_col(show.legend = FALSE) +
      geom_text(aes(label = label), vjust = -0.5, size = 3.5) +
      labs(x = "", y = "Percentage") +
      ggtitle(var) +
      ylim(0, y_limit) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, hjust = 0.5, margin = margin(b = 10))
      )
  })
  
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  
  if (length(plot_list) > 0) {
    pdf_file <- paste0(gene_name, "_", rs_id, "_GEpatientMetadata_BCIvariants_barplots.pdf")
    pdf(pdf_file, width = 10, height = 8)
    
    # Combinar plots y poner título general
    print(wrap_plots(plot_list, ncol = 2) +
      plot_annotation(
        title = paste0(rs_id, " (", length(df$participant_id), ")"),
        theme = theme(
          plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 20))
        )
      ))
    
    dev.off()
  }
}