#!/usr/bin/env Rscript
# Script: plotQCsmallVar

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("script.R <gene_name> <clean_table> <prot_file> <exon_file>")
}

clean_table     <- args[1] # "C:\\Users\\qp241615\\OneDrive - Queen Mary, University of London\\Documents\\4. Projects\\1. DHX34\\results\\DDX41\\DDX41_small_variants.rds"
gene_name       <- args[2] # "DDX41"
prot_file       <- args[3] # "C:\\Users\\qp241615\\OneDrive - Queen Mary, University of London\\Documents\\4. Projects\\1. DHX34\\data\\reference\\Protein\\DDX41.gff"
exon_file       <- args[4] # "C:\\Users\\qp241615\\OneDrive - Queen Mary, University of London\\Documents\\4. Projects\\1. DHX34\\data\\reference\\Exon\\DDX41.tsv"


# Message validation files
cat("Validating input files for gene:", gene_name, "\n")
missing_files <- c()
if (!file.exists(clean_table)) missing_files <- c(missing_files, paste("Clean table:", clean_table))
if (!file.exists(prot_file)) missing_files <- c(missing_files, paste("Protein file:", prot_file))
if (!file.exists(exon_file)) missing_files <- c(missing_files, paste("Exon file:", exon_file))

if (length(missing_files) > 0) {
  cat("ERROR: Missing files detected:\n")
  cat(paste(missing_files, collapse = "\n"), "\n")
  stop("Cannot proceed with missing input files")
}

cat("All input files validated successfully\n")



# Libraries  -------------------------------------------------------------------
library(tidyverse)
library(rtracklayer)
library(trackViewer)
library(GenomicRanges)
library(paletteer)

# Settings ----------------------------------------------------------------------
set.seed(23)

# Palettete ---------------------------------------------------------------------
my_pal <- c(paletteer_d("RColorBrewer::Set1"),paletteer_d("RColorBrewer::Set3"))
pastel_colors <- c("#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", "#BAE1FF", "#D5BAFF", 
                   "#FFC8E1", "#C2F0FC","#C5FAD5", "#FFDAC1", "#E2F0CB", "#C3B1E1", "#F8B195", 
                   "#F67280","#F6E2B3", "#B5EAD7","#FFABAB", "#D0F4DE", "#A0CED9", "#FFCBF2")


# Load Data ---------------------------------------------------------------------
## Variants
variants_table <- readRDS(clean_table)

## Protein gff
protein_info <- rtracklayer::import(prot_file)

## Exon tsv
exon_info <- read_tsv(exon_file)


# Visual Plots -----------------------------------------------------------------
## Distribution of MAF_variants ------
pdf(paste0(gene_name, "_MAF_Distribution.pdf"))
ggplot(variants_table, aes(x = log10(MAF_variant + 1e-6))) +
      geom_histogram(bins = 50, fill = "black") +
          labs(title = paste0(gene_name, " MAF Distribution"), x = "MAF (log10) + 1e-6",
                     y = "Number of variants") +
      theme_bw()
dev.off()


## MAF vs IMPACT -------
pdf(paste0(gene_name,"_MAFvsIMPACT.pdf"), width=14, height=8)
ggplot(variants_table, aes(x = log10(MAF_variant + 1e-6), y = IMPACT_annotation, 
                  color = CLIN_SIG_annotation, shape = SYMBOL_annotation)) +
      geom_jitter(width = 0.1, height = 0.1, size = 3, alpha = 0.8) +
      scale_x_continuous(name = "log10(MAF + 1e-6)") +
      scale_y_discrete(name = "IMPACT") +
      scale_shape_manual(values = c(16, 17)) +  
      theme_minimal(base_size = 14) +
      theme(legend.position = "right",plot.title = element_text(face="bold", size=16)) +
      labs(title = "MAF vs. IMPACT", color = "ClinVar",shape = "Gene")
dev.off()


## IMPACT Frequency -------
pdf(paste0(gene_name,"_IMPACT_Frequency.pdf"))
print(ggplot(variants_table, aes(x = IMPACT_annotation)) + 
      geom_bar(fill = "#0072B2", color = "black", width = 0.7) +
      geom_text(stat="count", aes(label= after_stat(count)),vjust=-0.5, color="black", size=5) +
      theme_minimal(base_size = 14) +
      labs(title = paste0("Frequency of ", gene_name, " by IMPACT category"),
        x = "IMPACT Category", y = "Number of variants") +
      theme(plot.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(angle = 30, hjust = 1)))
dev.off()

## Density Variants -------
# Add a check to ensure there are at least two points in each group before calculating density
valid_groups <- variants_table %>% 
  group_by(IMPACT_annotation) %>% 
  filter(n() > 1)  # Filter groups with more than one point

if (nrow(valid_groups) > 0) {
  densities <- valid_groups %>% 
    group_by(IMPACT_annotation) %>% 
    group_map(~ density(.x$POS_variant)$y)

  # Filter out NA values from densities
  max_density <- max(unlist(densities), na.rm = TRUE)
  if (is.infinite(max_density)) {
    max_density <- 0  # Set to 0 if no valid densities are found
  }
} else {
  cat("Warning: Not enough data points for density calculation. Skipping.")
  max_density <- 0
}

pdf(paste0(gene_name,"_Density_variants.pdf"), width=14, height=8)
ggplot(variants_table, aes(x = POS_variant, color = IMPACT_annotation, fill = IMPACT_annotation)) +
    geom_density(bw=500,alpha = 0.3) +  # Bandwidth of 1kB
    geom_rect(data = exon_info,
              aes(xmin = ExonStart, xmax = ExonEnd, ymin = -max_density*0.03, ymax = max_density*0.03),
              fill = "#00441B", alpha = 1, inherit.aes = FALSE) +
    geom_text(data = exon_info,
              aes(x = (ExonStart + ExonEnd) / 2, y = -max_density*0.03 * 3, label = seq_len(nrow(exon_info))),
              size = 3.5, inherit.aes = FALSE) +
    scale_color_manual(values =paletteer_d("RColorBrewer::RdYlGn")[c(1,3,5,10)]) +
    scale_fill_manual(values = paletteer_d("RColorBrewer::RdYlGn")[c(1,3,5,10)]) +
    facet_wrap(~ CHROM_variant, scales = "free_x") +
    labs(title = paste0("Density distribution of ",gene_name ," variants in genomic positions by impact annotation"),
        x = "Genomic Position (exons in dark green)", y = "Density  bw = 0.5kb",color = "Impact", fill = "Impact") +
    theme_minimal()
dev.off()

    
# Lollipop Plot Function -----------------------------------------------------
# Function to create trackViewer lollipop plot
create_lollipop2_plot <- function(df_data, gene_name, subset_label, protein_info) {
  
  # Prepare features
  feature2plot <- c("Domain", "Motif", "Region", "Zinc finger")
  features <- protein_info[protein_info$type %in% feature2plot, ]
  
  if(length(features) == 0) {
    cat("No protein features available for Lollipop plot\n")
    return()
  }
  
  # Prepare colors and labels for features
  features$Note <- trimws(features$Note)
  note_colors <- setNames(my_pal[1:length(unique(features$Note))], unique(features$Note))
  features$fill <- note_colors[features$Note]
  names(features) <- features$Note
  
  # Create variants GRanges
  sample.gr <- GRanges(seqnames = paste0(unique(seqnames(features))), 
                       ranges = IRanges(start = df_data$Protein_pos_start, width = 1, 
                                       names = df_data$LabelVarPlot))
  
  # Prepare variant properties
  sample.gr$score <- df_data$NS_variant
  sample.gr$Consequence_annotation <- df_data$Consequence_annotation
  sample.gr$color <- setNames(pastel_colors, unique(sample.gr$Consequence_annotation))[sample.gr$Consequence_annotation]
  legends <- list(labels=unique(sample.gr$Consequence_annotation), fill=unique(sample.gr$color))
  sample.gr$shape <- ifelse(grepl("^rs", names(sample.gr)), "circle", "diamond")
  
  # Prepare for plotting
  sample.gr.rot <- sample.gr
  sample.gr.rot$label.parameter.rot <- 45
  
  # Create filename and plot
  filename <- paste0(gene_name, "_Lollipop_", subset_label, ".pdf")

  pdf(filename, width=10, height=5)
  lolliplot(sample.gr.rot, features, legend=legends, ylab="Num. of Samples",
            yaxis.gp = gpar(fontsize=15), xaxis.gp = gpar(fontsize=15))
  grid.text(paste0(subset_label, " Variants of ", gene_name), 
                   x=.5, y=.95, gp=gpar(cex=1.5, fontface="bold"))
  dev.off()
  
  cat("Generated:", filename, "\n")
}

# Configuration for variant subsets -------------------------------------------
# Define which combinations to analyze
variant_subsets <- list(
  # High priority combinations
  list(
    name = "HIGH_All",
    filters = list(IMPACT_annotation = "HIGH"),
    description = "High_Impact_All"
  ),
  list(
    name = "HIGH_Canonical",
    filters = list(IMPACT_annotation = "HIGH", CANONICAL_annotation = "YES"),
    description = "High_Impact_Canonical"
  ),
  list(
    name = "HIGH_NonCanonical", 
    filters = list(IMPACT_annotation = "HIGH"),
    special_filters = list(CANONICAL_annotation = "not_YES"),  # This will include NA and other values
    description = "High_Impact_NonCanonical"
  ),
  list(
    name = "MODERATE_All",
    filters = list(IMPACT_annotation = "MODERATE"),
    description = "Moderate_Impact_All"
  ),
  list(
    name = "MODERATE_Canonical",
    filters = list(IMPACT_annotation = "MODERATE", CANONICAL_annotation = "YES"), 
    description = "Moderate_Impact_Canonical"
  )#,
  # Lower priority combinations (uncomment if needed)
  # list(
  #   name = "MODERATE_NonCanonical",
  #   filters = list(IMPACT_annotation = "MODERATE"),
  #   special_filters = list(CANONICAL_annotation = "not_YES"),
  #   description = "Moderate_Impact_NonCanonical"
  # ),
  # list(
  #   name = "LOW_All", 
  #   filters = list(IMPACT_annotation = "LOW"),
  #   description = "Low_Impact_All"
  # )
)

# Generate Lollipop Plots ------------------------------------------------------
cat("Generating lollipop plots for prioritized variant subsets...\n")

for (subset_config in variant_subsets) {
  
  # Apply standard filters
  df_subset <- variants_table
  for (column in names(subset_config$filters)) {
    if (column %in% colnames(variants_table)) {
      df_subset <- df_subset %>% 
        filter(!!sym(column) == subset_config$filters[[column]])
    } else {
      cat("Warning: Column", column, "not found in data. Skipping subset", subset_config$name, "\n")
      next
    }
  }
  
  # Apply special filters (for more complex filtering like "not YES")
  if (!is.null(subset_config$special_filters)) {
    for (column in names(subset_config$special_filters)) {
      if (column %in% colnames(variants_table)) {
        filter_value <- subset_config$special_filters[[column]]
        if (filter_value == "not_YES") {
          # Include everything except "YES" (includes NA, "NO", etc.)
          df_subset <- df_subset %>% 
            filter(is.na(!!sym(column)) | !!sym(column) != "YES")
        }
        # Add more special filter types here if needed
      } else {
        cat("Warning: Column", column, "not found in data. Skipping special filter\n")
      }
    }
  }
  
  # Additional base filters for lollipop plots
  df_subset <- df_subset %>% 
    filter(!is.na(Protein_pos_start), !is.na(NS_variant), NS_variant > 0)
  
  # Check if we have data to plot
  if (nrow(df_subset) == 0) {
    cat("No variants found for subset:", subset_config$description, "\n")
    next
  }
  

  # Generate lollipop plot
  cat("Processing subset:", subset_config$description, "(", nrow(df_subset), "variants )\n")
  create_lollipop2_plot(df_subset, gene_name, subset_config$description, protein_info)
}

cat("Plots generated successfully for gene:", gene_name, "\n")
