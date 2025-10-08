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
library(grid)
library(gridExtra)


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
metadata <- readRDS(metadata_file)


# GE samples with BCI variants --------------------------------------------------
ge_samples_with_bci_var <- comparison_data %>% 
  select(Variant_Name, Het_samples, Hom_samples, Hemi_samples) %>% 
  mutate(rs_id = str_extract(Variant_Name, "^[^,]+")) %>%
  filter(!is.na(rs_id)) %>%
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
  samples_vector <- ge_samples_with_bci_var[[rs_id]]
    write.table(samples_vector, 
                file = paste0(gene_name, "_GEsamples_with_BCIvariants_", rs_id, ".csv"),
                sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
}



# GE patients metadata with BCI variants ----------------------------------------
ge_patients_with_bci_var <- lapply(ge_samples_with_bci_var, function(samples) {
  metadata %>% filter(plate_key %in% samples)
})

# Save in csv files
for (rs_id in names(ge_patients_with_bci_var)) {
  df <- ge_patients_with_bci_var[[rs_id]]
    write.table(df, 
                file = paste0(gene_name, "_GEpatientMetadata_BCIvariants_", rs_id, ".csv"),
                sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
    
    df2 <- df %>% select(participant_id)
    write.table(df2, 
                file = paste0(gene_name, "_ParticipantsID_", rs_id, ".csv"),
                sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
}


# Metadata Distribution Plots ------------------------------------------------
# Function to create metadata plots for a given dataset
create_metadata_plots <- function(variant_metadata, rs_id, gene_name) {
  if (nrow(variant_metadata) == 0) {
    cat("No data available for", rs_id, "\n")
    return(NULL)
  }
  
  cat("Processing metadata for", nrow(variant_metadata), "participants with", rs_id, "variant...\n")
  
  # Ensure unique participants and add diagnosis age
  variant_metadata <- variant_metadata %>%
    distinct(participant_id, .keep_all = TRUE) %>%
    mutate(diagnosis_age = coalesce(
      as.numeric(as.character(cancer_diagnosis_age)),
      as.numeric(as.character(rare_disease_diagnosis_age))
    ))
  
  # Variables of interest for barplots
  vars_of_interest <- c("participant_type", "affection_status", "programme",
                        "yob", "diagnosis_age", "participant_karyotyped_sex",  
                        "genetically_inferred_ancestry_thr", "normalised_disease_group", 
                        "normalised_specific_disease")
  
  # Create barplots for each variable
  plot_list <- list()
  
  for (var in vars_of_interest) {
    if (!var %in% colnames(variant_metadata)) {
      cat("Warning: Variable", var, "not found in metadata. Skipping...\n")
      next
    }
    
    # Handle variables depending on type
    if (is.numeric(variant_metadata[[var]]) || var %in% c("yob", "diagnosis_age")) {
      # Create bins for numeric variables
      df_count <- variant_metadata %>%
        filter(!is.na(.data[[var]]), .data[[var]] != "") %>%
        mutate(
          numeric_var = as.numeric(as.character(.data[[var]])),
          var_binned = case_when(
            var == "yob" ~ paste0(floor(numeric_var/10)*10, "s"),
            var == "diagnosis_age" ~ paste0(floor(numeric_var/10)*10, "s"),
            TRUE ~ as.character(numeric_var)
          ),
          order_value = case_when(
            var == "yob" ~ floor(numeric_var/10)*10,
            var == "diagnosis_age" ~ floor(numeric_var/10)*10,
            TRUE ~ numeric_var
          )
        ) %>%
        filter(!is.na(numeric_var)) %>%
        count(var_binned, order_value) %>%
        mutate(
          pct = n / sum(n) * 100,
          label = ifelse(pct < 5, "<5%", paste0(round(pct), "%"))
        ) %>%
        arrange(order_value)
    } else {
      # Handle categorical variables
      df_count <- variant_metadata %>%
        filter(!is.na(.data[[var]]), .data[[var]] != "") %>%
        count(.data[[var]]) %>%
        mutate(
          pct = n / sum(n) * 100,
          var_binned = .data[[var]],
          order_value = NA_real_
        ) %>%
        arrange(desc(n))
      
      # For disease group variables, keep only top 10 and combine rest
      if (var %in% c("normalised_disease_group", "normalised_specific_disease")) {
        if (nrow(df_count) > 10) {
          top_10 <- df_count[1:10, ]
          rest_count <- sum(df_count[11:nrow(df_count), ]$n)
          rest_pct <- sum(df_count[11:nrow(df_count), ]$pct)
          
          rest_row <- data.frame(
            tmp = "REST, other groups or combinations",
            n = rest_count,
            pct = rest_pct,
            var_binned = "REST, other groups or combinations",
            order_value = NA_real_,
            stringsAsFactors = FALSE
          )
          names(rest_row)[1] <- var
          
          df_count <- rbind(top_10, rest_row)
        }
      }
      
      df_count <- df_count %>%
        mutate(label = ifelse(pct < 5, "<5%", paste0(round(pct), "%")))
    }
    
    if (nrow(df_count) == 0) next
    
    # Truncate labels that are 50+ characters
    df_count <- df_count %>%
      mutate(
        var_binned = ifelse(
          nchar(as.character(var_binned)) >= 50,
          paste0(substr(as.character(var_binned), 1, 47), "..."),
          as.character(var_binned)
        )
      )
    
    max_label_length <- max(nchar(as.character(df_count$var_binned)), na.rm = TRUE)
    max_pct <- max(df_count$pct)
    y_limit <- max_pct * 1.15
    
    # Create plot with appropriate ordering
    if (var %in% c("yob", "diagnosis_age")) {
      p <- ggplot(df_count, aes(x = reorder(var_binned, order_value), y = pct, fill = var_binned)) +
        geom_col(show.legend = FALSE, alpha = 0.8) +
        geom_text(aes(label = label), vjust = -0.5, size = 3) +
        labs(x = "", y = "Perc. of Participants") +
        ggtitle(paste(var, "(n =", sum(df_count$n), ")")) +
        ylim(0, y_limit) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          plot.title = element_text(size = 11, hjust = 0.5, margin = margin(b = 10)),
          panel.grid.minor = element_blank()
        )
    } else if (var %in% c("normalised_disease_group", "normalised_specific_disease")) {
      df_count <- df_count %>%
        mutate(
          is_rest = var_binned == "REST, other groups or combinations",
          sort_order = ifelse(is_rest, Inf, -pct)
        ) %>%
        arrange(sort_order)
      
      p <- ggplot(df_count, aes(x = factor(var_binned, levels = var_binned), y = pct, fill = var_binned)) +
        geom_col(show.legend = FALSE, alpha = 0.8) +
        geom_text(aes(label = label), vjust = -0.5, size = 3) +
        labs(x = "", y = "Perc. of Participants") +
        ggtitle(paste(var, "(n =", sum(df_count$n), ")")) +
        ylim(0, y_limit) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
          plot.title = element_text(size = 11, hjust = 0.5, margin = margin(b = 10)),
          panel.grid.minor = element_blank()
        )
    } else {
      p <- ggplot(df_count, aes(x = reorder(var_binned, -pct), y = pct, fill = var_binned)) +
        geom_col(show.legend = FALSE, alpha = 0.8) +
        geom_text(aes(label = label), vjust = -0.5, size = 3) +
        labs(x = "", y = "Perc. of Participants") +
        ggtitle(paste(var, "(n =", sum(df_count$n), ")")) +
        ylim(0, y_limit) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          plot.title = element_text(size = 11, hjust = 0.5, margin = margin(b = 10)),
          panel.grid.minor = element_blank()
        )
    }
    
    plot_list[[var]] <- list(plot = p, max_label_length = max_label_length)
  }
  
  # Remove NULL plots
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  
  if (length(plot_list) > 0) {
    # Separate plots by label length
    short_label_plots <- list()
    long_label_plots <- list()
    
    for (var_name in names(plot_list)) {
      if (var_name %in% c("normalised_disease_group", "normalised_specific_disease") || plot_list[[var_name]]$max_label_length > 25) {
        long_label_plots[[var_name]] <- plot_list[[var_name]]$plot
      } else {
        short_label_plots[[var_name]] <- plot_list[[var_name]]$plot
      }
    }
    
    # Create PDF with all barplots
    pdf_file <- paste0(gene_name, "_GEpatientMetadata_BCIvariants_", rs_id, ".pdf")
    pdf(pdf_file, width = 12, height = 16)
    
    grid.newpage()
    
    # Add main title
    grid.text(paste0("Participant Distribution for ", gene_name, " - ", rs_id),
              x = 0.5, y = 0.97, 
              gp = gpar(fontsize = 16, fontface = "bold"))
    
    # FIRST GRID: Short label variables
    if (length(short_label_plots) > 0) {
      grid.text(paste0("General Participant Characteristics (", nrow(variant_metadata), ")"),
                x = 0.5, y = 0.93,
                gp = gpar(fontsize = 14, fontface = "bold"))
      
      n_short_plots <- length(short_label_plots)
      ncol_short <- 3
      nrow_short <- ceiling(n_short_plots / ncol_short)
      
      plot_width_short <- 1 / ncol_short
      plot_height_short <- 0.48 / nrow_short 
      
      for (i in seq_along(short_label_plots)) {
        row_idx <- ceiling(i / ncol_short)
        col_idx <- ((i - 1) %% ncol_short) + 1
        
        x_pos <- (col_idx - 0.5) * plot_width_short
        y_pos <- 0.89 - (row_idx - 0.5) * plot_height_short 
        
        vp <- viewport(x = x_pos, y = y_pos, 
                       width = plot_width_short * 0.95, 
                       height = plot_height_short * 0.95)
        
        print(short_label_plots[[i]], vp = vp)
      }
    }
    
    # SECOND GRID: Long label variables
    if (length(long_label_plots) > 0) {
      grid.text("Disease Group Characteristics",
                x = 0.5, y = 0.38,
                gp = gpar(fontsize = 14, fontface = "bold"))
      
      n_long_plots <- length(long_label_plots)
      ncol_long <- 2
      nrow_long <- ceiling(n_long_plots / ncol_long)
      
      plot_width_long <- 1 / ncol_long
      plot_height_long <- 0.32 / nrow_long
      
      for (i in seq_along(long_label_plots)) {
        row_idx <- ceiling(i / ncol_long)
        col_idx <- ((i - 1) %% ncol_long) + 1
        
        x_pos <- (col_idx - 0.5) * plot_width_long
        y_pos <- 0.34 - (row_idx - 0.5) * plot_height_long
        
        vp <- viewport(x = x_pos, y = y_pos, 
                       width = plot_width_long * 0.95, 
                       height = plot_height_long * 0.95)
        
        print(long_label_plots[[i]], vp = vp)
      }
    }
    
    dev.off()
    cat("Generated:", pdf_file, "\n")
  }
}

# Process each rs_id separately
for (rs_id in names(ge_patients_with_bci_var)) {
  cat("\n--- Processing", rs_id, "variant ---\n")
  df <- ge_patients_with_bci_var[[rs_id]]
  create_metadata_plots(df, rs_id, gene_name)
}

