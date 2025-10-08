#!/usr/bin/env Rscript
# Script: plotFilteredStructVar

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("script.R <filtered_table> <gene_name> <prot_file> <exon_file> <part_metadata_file>")
}

filtered_table  <- args[1] # "C:\\Users\\qp241615\\OneDrive - Queen Mary, University of London\\Documents\\4. Projects\\1. DHX34\\results\\DHX34\\DHX34_small_variants.rds"
gene_name       <- args[2] # "DHX34"
prot_file       <- args[3] # "C:\\Users\\qp241615\\OneDrive - Queen Mary, University of London\\Documents\\4. Projects\\1. DHX34\\data\\reference\\Protein\\DHX34.gff"
exon_file       <- args[4] # "C:\\Users\\qp241615\\OneDrive - Queen Mary, University of London\\Documents\\4. Projects\\1. DHX34\\data\\reference\\Exon\\DHX34.tsv"
p_metadata_file <- args[5] # "C:\\Users\\qp241615\\OneDrive - Queen Mary, University of London\\Documents\\4. Projects\\1. DHX34\\results\\DHX34\\DHX34_small_variants_participantMetadata.tsv"

# Message validation files
cat("Validating input files for gene:", gene_name, "\n")
missing_files <- c()
if (!file.exists(filtered_table)) missing_files <- c(missing_files, paste("Filtered table:", filtered_table))
if (!file.exists(prot_file)) missing_files <- c(missing_files, paste("Protein file:", prot_file))
if (!file.exists(exon_file)) missing_files <- c(missing_files, paste("Exon file:", exon_file))
if (!file.exists(p_metadata_file)) missing_files <- c(missing_files, paste("Participant metadata file:", p_metadata_file))

if (length(missing_files) > 0) {
  cat("ERROR: Missing files detected:\n")
  cat(paste(missing_files, collapse = "\n"), "\n")
  stop("Cannot proceed with missing input files")
}

cat("All input files validated successfully\n")


# Libraries  -------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(circlize)
library(paletteer)
library(grid)


# Settings ----------------------------------------------------------------------
set.seed(23)


# Load Data ---------------------------------------------------------------------
## Structural Variants
variants_table <- readRDS(filtered_table)

## Protein gff
protein_info <- rtracklayer::import(prot_file)

## Exon tsv
exon_info <- read_tsv(exon_file)

## Participant Metadata
metadata_info <- readRDS(p_metadata_file)



# Plots Functions ---------------------------------------------------------
generate_barplot <- function(data, gene, variant_orig) {
  pdf(paste0(gene, "_Filtered_", variant_orig, "_Frec_SVtype.pdf"), width = 10, height = 8)
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
  pdf(paste0(gene, "_Filtered_", variant_orig, "_Circle_Translocation.pdf"), width = 8, height = 8)
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
  pdf(paste0(gene, "_Filtered_", variant_orig, "_CNV_Size_CN.pdf"), width = 10, height = 8)
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
  pdf(paste0(gene, "_Filtered_", variant_orig, "_Inversions_Segments.pdf"), width = 14, height = 10)
  
  # Get gene boundaries from exon info
  gene_start <- min(exon_info$ExonStart, na.rm = TRUE)
  gene_end <- max(exon_info$ExonEnd, na.rm = TRUE)
  
  plot_data <- data %>%
    mutate(start_pos = as.numeric(POS), end_pos = as.numeric(INFO_END), quality = as.numeric(QUAL_SVonly), 
          unique_inversion = paste0(CHROM, "_", POS, "_", INFO_END)) %>%
    group_by(unique_inversion, start_pos, end_pos) %>%
    summarise(count = n(), avg_quality = mean(quality, na.rm = TRUE), .groups = "drop") %>%
    mutate(inversion_length = end_pos - start_pos,
          inversion_label = paste0(start_pos, "_", end_pos, "_n=", count)) %>%
    arrange(desc(inversion_length)) %>%  # Sort by length (longest first)
    mutate(variant_id = row_number())    # Assign IDs based on length order
  
  # PAGE 1: Full view
  p1 <- ggplot(plot_data, aes(y = variant_id)) +
    # Add gene bar
    geom_segment(aes(x = gene_start, xend = gene_end, y = 0, yend = 0), 
                 color = "black", size = 3, alpha = 0.8) +
    geom_text(aes(x = (gene_start + gene_end)/2, y = 0, label = paste(gene, "Gene")), 
              vjust = 2, size = 4) +
    # Add inversion segments
    geom_segment(aes(x = start_pos, xend = end_pos, yend = variant_id, color = avg_quality),
                size = 1.2, alpha = 0.8) +
    geom_text(aes(x = (start_pos + end_pos)/2, label = inversion_label), 
              vjust = -0.5, size = 2.5, angle = 0) +
    scale_color_gradient(low = "#0066CC", high = "#FF0000", name = "Avg Quality") +
    theme_minimal() +
    labs(title = paste0("Inversions in ", variant_orig, " variants for ", gene),
      x = "Genomic Position (bp)", y = "Unique Inversions") +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
      plot.title = element_text(size = 16), legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)) +
    scale_x_continuous(labels = scales::comma_format(), 
                       limits = c(min(gene_start, min(plot_data$start_pos)), 
                                 max(gene_end, max(plot_data$end_pos)))) +
    scale_y_continuous(breaks = NULL, limits = c(-0.5, max(plot_data$variant_id) + 1))
  
  print(p1)
  
  # PAGE 2: Zoomed view (±1MB around gene)
  zoom_start <- gene_start - 1000000
  zoom_end <- gene_end + 1000000
  
  plot_data_zoom <- plot_data %>%
    filter(start_pos <= zoom_end & end_pos >= zoom_start) %>%
    mutate(start_pos_clipped = pmax(start_pos, zoom_start),
          end_pos_clipped = pmin(end_pos, zoom_end),
      inversion_label_zoom = case_when( # Update labels to show if clipped
        start_pos < zoom_start & end_pos > zoom_end ~ paste0("...", start_pos, "_", end_pos, "_n=", count, "..."),
        start_pos < zoom_start ~ paste0("...", start_pos, "_", end_pos, "_n=", count),
        end_pos > zoom_end ~ paste0(start_pos, "_", end_pos, "_n=", count, "..."),
        TRUE ~ inversion_label
      )
    )
  
  if (nrow(plot_data_zoom) > 0) {
    p2 <- ggplot(plot_data_zoom, aes(y = variant_id)) +
      # Add gene bar at the bottom
      geom_segment(aes(x = gene_start, xend = gene_end, y = 0, yend = 0), 
                   color = "black", size = 3, alpha = 0.8) +
      geom_text(aes(x = (gene_start + gene_end)/2, y = 0, label = paste(gene, "Gene")), 
                vjust = 2, size = 4) +
      # Add inversion segments (using clipped coordinates)
      geom_segment(aes(x = start_pos_clipped, xend = end_pos_clipped, yend = variant_id, color = avg_quality),
                  size = 1.5, alpha = 0.8) +
      geom_text(aes(x = (start_pos_clipped + end_pos_clipped)/2, label = inversion_label_zoom), 
                vjust = -0.5, size = 3, angle = 0) +
      scale_color_gradient(low = "#0066CC", high = "#FF0000", name = "Avg Quality") +
      theme_minimal() +
      labs(title = paste0("Inversions in ", variant_orig, " variants for ", gene, " - Gene Region (±1MB)"),
        x = "Genomic Position (bp)", y = "Unique Inversions") +
      theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 16), legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
      scale_x_continuous(labels = scales::comma_format(), 
                         limits = c(zoom_start, zoom_end)) +
      scale_y_continuous(breaks = NULL, limits = c(-0.5, max(plot_data_zoom$variant_id) + 1))
    
    print(p2)
  }   
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
# Function to create metadata plots for a given dataset
create_metadata_plots <- function(variant_data, title_suffix, filename_suffix, variant_orig) {
  if (nrow(variant_data) == 0) {
    cat("No data available for", title_suffix, "\n")
    return(NULL)
  }
  
  cat("Processing metadata for", nrow(variant_data), title_suffix, "variants...\n")
  
  # Extract all unique samples from the SAMPLE column
  all_samples <- variant_data %>%
    pull(SAMPLE) %>%
    unique() %>%
    .[!is.na(.) & . != ""]
  
  cat("Found", length(all_samples), "unique samples with", title_suffix, "variants\n")
  
  # Match samples to participant metadata and ensure unique participants
  variant_metadata <- metadata_info %>%
    filter(plate_key %in% all_samples) %>%
    distinct(participant_id, .keep_all = TRUE) %>%  # Keep only unique participants
    mutate(diagnosis_age = coalesce(
        as.numeric(as.character(cancer_diagnosis_age)),
        as.numeric(as.character(rare_disease_diagnosis_age))
      )
    )
  
  cat("Matched to", nrow(variant_metadata), "unique participants\n")
  
  if (nrow(variant_metadata) == 0) {
    cat("No participant metadata found for", title_suffix, "variants\n")
    return(NULL)
  }
  
  # Variables of interest for barplots
  vars_of_interest <- c("participant_type", "affection_status", "programme",
                        "yob", "diagnosis_age", "participant_karyotyped_sex",  
                        "genetically_inferred_ancestry_thr","normalised_disease_group", "normalised_disease_sub_group")
  
  # Create barplots for each variable
  plot_list <- list()
  
  for (var in vars_of_interest) {
    if (!var %in% colnames(variant_metadata)) {
      cat("Warning: Variable", var, "not found in metadata. Skipping...\n")
      next
    }
    
    # Handle variables depending on type
    if (is.numeric(variant_metadata[[var]]) || var %in% c("yob", "diagnosis_age")) {
      # Create bins for numeric variables - ensure proper numeric conversion
      df_count <- variant_metadata %>%
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
      # Handle categorical variables
      df_count <- variant_metadata %>%
        filter(!is.na(.data[[var]]), .data[[var]] != "") %>%
        count(.data[[var]]) %>%
        mutate(
          pct = n / sum(n) * 100,
          var_binned = .data[[var]],
          order_value = NA_real_  # No ordering for categorical
        ) %>%
        arrange(desc(n))
      
      # For disease group variables, keep only top 10 and combine rest
      if (var %in% c("normalised_disease_group", "normalised_disease_sub_group")) {
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
      
      # Add percentage labels
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
    # Separate plots by label length
    short_label_plots <- list()
    long_label_plots <- list()
    
    for (var_name in names(plot_list)) {  
      if (plot_list[[var_name]]$max_label_length > 25) {
        long_label_plots[[var_name]] <- plot_list[[var_name]]$plot
      } else {
        short_label_plots[[var_name]] <- plot_list[[var_name]]$plot
      }
    }
    
    # Create PDF with all barplots on the same page using grid functions
    pdf_file <- paste0(gene_name, "_Filtered_",variant_orig,"_PartMetadata_Barplots_", filename_suffix, ".pdf")
    pdf(pdf_file, width = 12, height = 16)
    
    # Create a single page with both grids
    grid.newpage()
    
    # Add main title for the entire page
    grid.text(paste0("Participant Distribution for ", gene_name, " (", nrow(variant_data), " ", title_suffix, " Filtered Variants)"),
              x = 0.5, y = 0.97, 
              gp = gpar(fontsize = 16, fontface = "bold"))
    
    # FIRST GRID: Short label variables (≤25 characters) - Upper part
    if (length(short_label_plots) > 0) {
      # Add section title for short labels - moved up
      grid.text(paste0("General Participant Characteristics (", nrow(variant_metadata), ")"),
                x = 0.5, y = 0.93,  # Moved from 0.89 to 0.93
                gp = gpar(fontsize = 14, fontface = "bold"))
      
      # Calculate grid layout (3 columns for short labels)
      n_short_plots <- length(short_label_plots)
      ncol_short <- 3
      nrow_short <- ceiling(n_short_plots / ncol_short)
      
      # Define viewport dimensions for upper half - increased height
      plot_width_short <- 1 / ncol_short
      plot_height_short <- 0.48 / nrow_short 
      
      # Create viewports and print plots for short labels
      for (i in seq_along(short_label_plots)) {
        row_idx <- ceiling(i / ncol_short)
        col_idx <- ((i - 1) %% ncol_short) + 1
        
        # Calculate viewport position (upper half) - adjusted starting position
        x_pos <- (col_idx - 0.5) * plot_width_short
        y_pos <- 0.89 - (row_idx - 0.5) * plot_height_short 
        
        # Create viewport
        vp <- viewport(x = x_pos, y = y_pos, 
                       width = plot_width_short * 0.95, 
                       height = plot_height_short * 0.95)
        
        # Print plot in viewport
        print(short_label_plots[[i]], vp = vp)
      }
    }
    
    # SECOND GRID: Long label variables (>25 characters) - Lower part
    if (length(long_label_plots) > 0) {
      # Add section title for long labels - moved down
      grid.text("Disease Group Characteristics",
                x = 0.5, y = 0.38,  # Moved from 0.45 to 0.38
                gp = gpar(fontsize = 14, fontface = "bold"))
      
      # Calculate grid layout (2 columns for long labels to give more space)
      n_long_plots <- length(long_label_plots)
      ncol_long <- 2
      nrow_long <- ceiling(n_long_plots / ncol_long)
      
      # Define viewport dimensions for lower half - reduced height
      plot_width_long <- 1 / ncol_long
      plot_height_long <- 0.32 / nrow_long  # Reduced from 0.38 to 0.32
      
      # Create viewports and print plots for long labels
      for (i in seq_along(long_label_plots)) {
        row_idx <- ceiling(i / ncol_long)
        col_idx <- ((i - 1) %% ncol_long) + 1
        
        # Calculate viewport position (lower half) - moved down
        x_pos <- (col_idx - 0.5) * plot_width_long
        y_pos <- 0.34 - (row_idx - 0.5) * plot_height_long  # Moved from 0.41 to 0.34
        
        # Create viewport
        vp <- viewport(x = x_pos, y = y_pos, 
                       width = plot_width_long * 0.95, 
                       height = plot_height_long * 0.95)
        
        # Print plot in viewport
        print(long_label_plots[[i]], vp = vp)
      }
    }
    
    dev.off()
    
    cat("Generated:", pdf_file, "\n")
    cat("Short label variables (≤25 chars):", length(short_label_plots), "\n")
    cat("Long label variables (>25 chars):", length(long_label_plots), "\n")
  } else {
    cat("No valid plots generated for", title_suffix, "metadata\n")
  }
  
  return(variant_metadata)
}

for (variant_orig in unique(variants_table$variant_origin)) {
  # Subset the data for variant origin
  variants_table_VarOrg <- variants_table[variants_table$variant_origin == variant_orig, ]
  if (nrow(variants_table_VarOrg) == 0) {
    cat("No data available for variant origin:", variant_orig, "\n")
    next
  }

  # Get unique structural variant types
  sv_types <- unique(variants_table_VarOrg$INFO_SVTYPE)
  sv_types <- sv_types[!is.na(sv_types) & sv_types != ""]

  cat("Found structural variant types:", paste(sv_types, collapse = ", "), "\n")

  # Process each structural variant type separately
  for (sv_type in sv_types) {
    cat("\n--- Processing", sv_type, "variants ---\n")
    
    # Filter variants for this specific type
    sv_type_variants <- variants_table_VarOrg %>% 
      filter(INFO_SVTYPE == sv_type)
    
    # Create metadata plots for this SV type
    create_metadata_plots(
      variant_data = sv_type_variants,
      title_suffix = sv_type,
      filename_suffix = sv_type,
      variant_orig = variant_orig
    )
  }

  # Create combined analysis for all structural variants
  cat("\n--- Processing All Structural Variants ---\n")
  create_metadata_plots(
    variant_data = variants_table_VarOrg,
    title_suffix = "All Structural",
    filename_suffix = "AllStructuralVar",
    variant_orig = variant_orig
  )
}

cat("Plots generated successfully for gene:", gene_name, "\n")



