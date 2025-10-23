#!/usr/bin/env Rscript
# Script: plotQCstructVar.R

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("script.R <struct_variants> <gene_name> <part_metadata_file>")
}

struct_variants     <- args[1] # r"(C:\Users\qp241615\OneDrive - Queen Mary, University of London\Documents\4. Projects\1. DHX34\results\DDX41\DDX41_structural_variants.rds)"
gene_name           <- args[2] # "DDX41"
part_metadata_file  <- args[3] # "C:\\Users\\qp241615\\OneDrive - Queen Mary, University of London\\Documents\\4. Projects\\1. DHX34\\results\\DHX34\\DHX34_small_variants_participantMetadata.tsv"


# Message validation files
cat("Validating input files for gene:", gene_name, "\n")
missing_files <- c()
if (!file.exists(struct_variants)) missing_files <- c(missing_files, paste("Structural variants:", struct_variants))
if (!file.exists(part_metadata_file)) missing_files <- c(missing_files, paste("Participant metadata file:", part_metadata_file))

if (length(missing_files) > 0) {
  cat("ERROR: Missing files detected:\n", paste(missing_files, collapse = "\n"), "\n")
  stop("Validation failed due to missing files.")
}
cat("All input files validated successfully\n")


# Libraries  -------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(circlize)
library(paletteer)


# Settings ----------------------------------------------------------------------
set.seed(23)


# Load Data ---------------------------------------------------------------------
## Structural Variants
variants_table <- readRDS(struct_variants)

## Participant Metadata
metadata_info <- readRDS(part_metadata_file)


# Plots Functions ---------------------------------------------------------
generate_barplot <- function(data, gene, variant_orig) {
  pdf(paste0(gene, "_QC_", variant_orig, "_Frec_SVtype.pdf"), width = 10, height = 8)
  print(
    ggplot(data %>% count(INFO_SVTYPE), aes(x = INFO_SVTYPE, y = n)) +
      geom_col(fill = "#377EB8") +
      geom_text(aes(label = n), vjust = -0.5) +
      theme_minimal() +
      labs(title = paste0("Frequency of ", variant_orig, " variants for ", gene),
           x = "Type of Structural Variant", y = "Number of Structural Variants") +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
            axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
            plot.title = element_text(size = 16))
  )
  dev.off()
}

generate_circos_plot <- function(data, gene, variant_orig) {
  pdf(paste0(gene, "_QC_", variant_orig, "_Circle_Translocation.pdf"), width = 8, height = 8)
  circos.clear()
  circos.initializeWithIdeogram(species = "hg38")
  title(paste0("Translocation of ", gene, " variants in ", variant_orig, " samples"))
  legend("bottomleft", pch = 1, legend = paste0(nrow(data), " Translocation Events"))
  
  for (i in seq_len(nrow(data))) {
    circos.link(
      sector.index1 = data$chr1[i], point1 = data$pos1[i],
      sector.index2 = data$chr2[i], point2 = data$pos2[i],
      col = "#FF000080", border = NA
    )
  }
  dev.off()
}

generate_cnv_plot <- function(data, gene, variant_orig) {
  pdf(paste0(gene, "_QC_", variant_orig, "_CNV_Size_CN.pdf"), width = 10, height = 8)
  print(
    ggplot(data, aes(x = log10(SV_size), y = CN_CNVonly, color = CN_status)) +
      geom_point(alpha = 0.7, size = 3) +
      theme_minimal() +
      labs(
        title = paste0("Size vs Copy Number of ", variant_orig, " CNVs of ", gene),
        x = "log10(Size in bp)", y = "Copy Number"
      ) +
      scale_color_brewer(palette = "Set2", na.value = "grey50") +
      scale_y_continuous(breaks = seq(floor(min(data$CN_CNVonly, na.rm = TRUE)), 
                                      ceiling(max(data$CN_CNVonly, na.rm = TRUE)), by = 1))
  )
  dev.off()
}

generate_inversion_plot <- function(data, gene, variant_orig) {
  pdf(paste0(gene, "_QC_", variant_orig, "_Inversions_Segments.pdf"), width = 14, height = 10)
  
  plot_data <- data %>%
    mutate(start_pos = as.numeric(POS), end_pos = as.numeric(INFO_END), 
           quality = as.numeric(QUAL_SVonly),  # Asegurar que 'quality' es numérico
           unique_inversion = paste0(CHROM, "_", POS, "_", INFO_END)) %>%
    group_by(unique_inversion, start_pos, end_pos) %>%
    summarise(count = n(), avg_quality = mean(quality, na.rm = TRUE), .groups = "drop") %>%
    mutate(inversion_length = end_pos - start_pos) %>%
    arrange(desc(inversion_length)) %>%  # Sort by length (longest first)
    mutate(variant_id = row_number())    # Assign IDs based on length order
  
  # Single plot view
  p1 <- ggplot(plot_data, aes(y = variant_id)) +
    geom_segment(aes(x = start_pos, xend = end_pos, yend = variant_id, color = avg_quality),
                linewidth = 1.2, alpha = 0.8) +
    scale_color_gradient(low = "#0066CC", high = "#FF0000", name = "Avg Quality") +
    theme_minimal() +
    labs(title = paste0("Inversions in ", variant_orig, " variants for ", gene),
         x = "Genomic Position (bp)", y = "Unique Inversions") +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
          plot.title = element_text(size = 16), legend.title = element_text(size = 12),
          legend.text = element_text(size = 10)) +
    scale_x_continuous(labels = scales::comma_format()) +
    scale_y_continuous(breaks = NULL)
  
  print(p1)
  dev.off()
}


# Plots Generation -------------------------------------------------------------------
for (variant_orig in unique(variants_table$variant_origin)) {
  # Subset the data for variant origin
  variants_table_VarOrg <- variants_table[variants_table$variant_origin == variant_orig, ]
  if (nrow(variants_table_VarOrg) == 0) {
    cat("No data available for variant origin:", variant_orig, "\n")
    next
  }

  # Generate barplot
  generate_barplot(variants_table_VarOrg, gene_name, variant_orig)

  # Generate circos plot for translocations
  BND_df <- variants_table_VarOrg %>%
    filter(INFO_SVTYPE == "BND") %>%
    mutate(chr1 = CHROM, chr2 = str_extract(ALT, "chr[0-9XYM]+"), pos1 = POS,
           pos2 = as.numeric(str_extract(ALT, "(?<=:)\\d+"))
    ) %>%
    select(chr1, pos1, chr2, pos2)

  if (nrow(BND_df) == 0) {
    cat("No translocation data available for variant origin:", variant_orig, "\n")
  } else {
    generate_circos_plot(BND_df, gene_name, variant_orig)
  }

  # Generate CNV plot
  variants_table_VarOrg_CNV <- variants_table_VarOrg %>%
    filter(INFO_SVTYPE == "CNV") %>%
    mutate(
      SV_size = as.numeric(INFO_END) - as.numeric(POS),
      CN_status = case_when(
        CN_CNVonly > 2 ~ "Gain", CN_CNVonly < 2 ~ "Loss",
        CN_CNVonly == 2 ~ "Neutral", TRUE ~ NA_character_
      )
    )

  if (nrow(variants_table_VarOrg_CNV) == 0) {
    cat("No CNV data available for variant origin:", variant_orig, "\n")
  } else {
    generate_cnv_plot(variants_table_VarOrg_CNV, gene_name, variant_orig)
  }

  # Generate inversion plot
  variants_table_VarOrg_INV <- variants_table_VarOrg %>%
    filter(INFO_SVTYPE == "INV")

  if (nrow(variants_table_VarOrg_INV) == 0) {
    cat("No inversion data available for variant origin:", variant_orig, "\n")
  } else {
    generate_inversion_plot(variants_table_VarOrg_INV, gene_name, variant_orig)
  }
}

# Metadata Distribution Plots ------------------------------------------------
# Function to generate metadata plots for structural variants by variant origin (germline/somatic)
generate_metadata_plots <- function(variants_data, gene, variant_orig, metadata_info) {
  cat("Processing metadata for", variant_orig, "structural variants:", nrow(variants_data), "variants...\n")
  
  # Extract all samples from the variants
  # For structural variants, samples are in SAMPLE column
  if (!"SAMPLE" %in% colnames(variants_data)) {
    cat("Warning: SAMPLE column not found in structural variants data. Skipping metadata plots.\n")
    return()
  }

  all_samples <- variants_data %>%
    pull(SAMPLE) %>%
    paste(collapse = ",") %>%
    strsplit(",") %>%
    unlist() %>%
    unique() %>%
    trimws() %>%
    .[. != "" & !is.na(.)]
  
  cat("Found", length(all_samples), "unique samples with", variant_orig, "structural variants\n")
  
  if (length(all_samples) == 0) {
    cat("No samples found for", variant_orig, "variants. Skipping metadata plots.\n")
    return()
  }
  
  # Match samples to participant metadata and ensure unique participants
  filtered_metadata <- metadata_info %>%
    filter(plate_key %in% all_samples) %>%
    distinct(participant_id, .keep_all = TRUE) %>%  # Keep only unique participants
    mutate(diagnosis_age = coalesce(
        as.numeric(as.character(cancer_diagnosis_age)),
        as.numeric(as.character(rare_disease_diagnosis_age))
      )
    )
  
  if (nrow(filtered_metadata) == 0) {
    cat("No participant metadata found for", variant_orig, "variants. Skipping plots.\n")
    return()
  }
  
  # Variables of interest for barplots
  vars_of_interest <- c("participant_type", "affection_status", "programme",
                        "yob", "diagnosis_age", "participant_karyotyped_sex",
                        "genetically_inferred_ancestry_thr", "normalised_disease_group", "normalised_disease_sub_group")
  
  # Create barplots for each variable
  plot_list <- list()
  
  for (var in vars_of_interest) {
    if (!var %in% colnames(filtered_metadata)) {
      cat("Warning: Variable", var, "not found in metadata. Skipping...\n")
      next
    }
    
    # Handle variables depending on type
    if (is.numeric(filtered_metadata[[var]]) || var %in% c("yob", "diagnosis_age")) {
      # Create bins for numeric variables - ensure proper numeric conversion
      df_count <- filtered_metadata %>%
        filter(!is.na(.data[[var]]), .data[[var]] != "") %>%
        mutate(
          numeric_var = as.numeric(as.character(.data[[var]])),  # Force numeric conversion
          var_binned = case_when(
            var == "yob" ~ paste0(floor(numeric_var/10)*10, "s"),
            var == "diagnosis_age" ~ paste0(floor(numeric_var/10)*10, "s"),
            TRUE ~ as.character(numeric_var)
          ),
          # Create numeric ordering column for temporal variables
          order_value = case_when(
            var == "yob" ~ floor(numeric_var/10)*10,
            var == "diagnosis_age" ~ floor(numeric_var/10)*10,
            TRUE ~ numeric_var
          )
        ) %>%
        filter(!is.na(numeric_var)) %>%  # Remove rows where conversion failed
        count(var_binned, order_value) %>%
        mutate(
          pct = n / sum(n) * 100,
          label = ifelse(pct < 5, "<5%", paste0(round(pct), "%"))
        ) %>%
        arrange(order_value)  # Order by the numeric value instead of count
    } else {
      # For categorical variables, count and calculate percentages
      df_count <- filtered_metadata %>%
        filter(!is.na(.data[[var]]), .data[[var]] != "") %>%
        count(.data[[var]]) %>%
        rename(var_binned = 1) %>%
        mutate(
          pct = n / sum(n) * 100,
          label = ifelse(pct < 5, "<5%", paste0(round(pct), "%"))
        )
    }
    
    # Skip if no data
    if (nrow(df_count) == 0) {
      cat("No data for variable", var, "in", variant_orig, "variants\n")
      next
    }
    
    # Calculate maximum label length for this variable (after truncation)
    max_label_length <- max(nchar(as.character(df_count$var_binned)), na.rm = TRUE)
    
    # Calculate y-axis limit to accommodate labels
    max_pct <- max(df_count$pct)
    y_limit <- max_pct * 1.15
    
    # Create plot with appropriate ordering
    if (var %in% c("yob", "diagnosis_age")) {
      # For temporal variables, order by numeric value (chronological order)
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
    } else if (var %in% c("normalised_disease_group", "normalised_disease_sub_group")) {
      # Special handling for disease group variables - REST group first, then by frequency
      df_count <- df_count %>% 
        mutate(
          is_rest = var_binned == "REST, other groups or combinations",
          sort_order = ifelse(is_rest, Inf, -pct)  # REST gets highest value, others by negative pct (desc)
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
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),  # Changed to vertical labels
          plot.title = element_text(size = 11, hjust = 0.5, margin = margin(b = 10)),
          panel.grid.minor = element_blank()
        )
    } else {
      # For other categorical variables, order by frequency (descending)
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
    
    # Store plot with its maximum label length
    plot_list[[var]] <- list(plot = p, max_label_length = max_label_length)
  }
  
  # Remove NULL plots
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  
  if (length(plot_list) > 0) {
    # Generate individual plots for each variable type (similar to structural variant types)
    for (var_name in names(plot_list)) {
      pdf_file <- paste0(gene, "_QC_", variant_orig, "_PartMetadata_Barplots_", toupper(var_name), ".pdf")
      pdf(pdf_file, width = 10, height = 8)
      print(plot_list[[var_name]]$plot)
      dev.off()
    }
    
    # Also create a combined "All" plot similar to the main structural variant plots
    # Separate plots by label length for the combined plot
    short_label_plots <- list()
    long_label_plots <- list()
    
    for (var_name in names(plot_list)) {
      # Force both disease group variables to always go to long label plots
      if (var_name %in% c("normalised_disease_group", "normalised_disease_sub_group") || plot_list[[var_name]]$max_label_length > 25) {
        long_label_plots[[var_name]] <- plot_list[[var_name]]$plot
      } else {
        short_label_plots[[var_name]] <- plot_list[[var_name]]$plot
      }
    }
    
    # Create PDF with all barplots on the same page using grid functions
    pdf_file <- paste0(gene, "_QC_", variant_orig, "_PartMetadata_Barplots_AllStructuralVar.pdf")
    pdf(pdf_file, width = 12, height = 16)
    
    # Load grid library
    library(grid)
    
    # Create a single page with both grids
    grid.newpage()
    
    # Add main title for the entire page
    grid.text(paste0("Participant Distribution for ", gene, " ", str_to_title(variant_orig), " Structural Variants (", nrow(variants_data), " variants)"),
              x = 0.5, y = 0.97, 
              gp = gpar(fontsize = 16, fontface = "bold"))
    
    # FIRST GRID: Short label variables (≤25 characters) - Upper part
    if (length(short_label_plots) > 0) {
      # Add section title for short labels
      grid.text(paste0("General Participant Characteristics (", nrow(filtered_metadata), ")"),
                x = 0.5, y = 0.93,
                gp = gpar(fontsize = 14, fontface = "bold"))
      
      # Calculate grid layout (3 columns for short labels)
      n_short_plots <- length(short_label_plots)
      n_short_cols <- min(3, n_short_plots)
      n_short_rows <- ceiling(n_short_plots / n_short_cols)
      
      # Define viewport for short label plots - upper portion
      pushViewport(viewport(x = 0.5, y = 0.75, width = 0.95, height = 0.35))
      
      # Create grid layout
      pushViewport(viewport(layout = grid.layout(n_short_rows, n_short_cols)))
      
      # Plot each short label plot
      for (i in seq_along(short_label_plots)) {
        row <- ceiling(i / n_short_cols)
        col <- ((i - 1) %% n_short_cols) + 1
        pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
        grid.draw(ggplotGrob(short_label_plots[[i]]))
        popViewport()
      }
      popViewport()
      popViewport()
    }
    
    # SECOND GRID: Long label variables (>25 characters) - Lower part
    if (length(long_label_plots) > 0) {
      # Add section title for long labels
      start_y <- ifelse(length(short_label_plots) > 0, 0.48, 0.88)
      grid.text("Disease Categories and Long Labels",
                x = 0.5, y = start_y,
                gp = gpar(fontsize = 14, fontface = "bold"))
      
      # Calculate grid layout (2 columns for long labels)
      n_long_plots <- length(long_label_plots)
      n_long_cols <- min(2, n_long_plots)
      n_long_rows <- ceiling(n_long_plots / n_long_cols)
      
      # Define viewport for long label plots - lower portion
      viewport_y <- ifelse(length(short_label_plots) > 0, 0.25, 0.65)
      viewport_height <- ifelse(length(short_label_plots) > 0, 0.35, 0.7)
      
      pushViewport(viewport(x = 0.5, y = viewport_y, width = 0.95, height = viewport_height))
      
      # Create grid layout
      pushViewport(viewport(layout = grid.layout(n_long_rows, n_long_cols)))
      
      # Plot each long label plot
      for (i in seq_along(long_label_plots)) {
        row <- ceiling(i / n_long_cols)
        col <- ((i - 1) %% n_long_cols) + 1
        pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
        grid.draw(ggplotGrob(long_label_plots[[i]]))
        popViewport()
      }
      popViewport()
      popViewport()
    }
    
    dev.off()
    cat("Generated metadata plots for", variant_orig, "structural variants\n")
  }
}

# Generate metadata plots for both germline and somatic variants
# Filter for germline variants
germline_variants <- variants_table %>% filter(variant_origin  == "germline")
if (nrow(germline_variants) > 0) {
  generate_metadata_plots(germline_variants, gene_name, "germline", metadata_info)
}

# Filter for somatic variants  
somatic_variants <- variants_table %>% filter(variant_origin  == "somatic")
if (nrow(somatic_variants) > 0) {
  generate_metadata_plots(somatic_variants, gene_name, "somatic", metadata_info)
}

cat("Plots generated successfully for gene:", gene_name, "\n")


