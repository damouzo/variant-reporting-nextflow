#!/usr/bin/env Rscript
# Script: plotQCstructVar.R

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("script.R <struct_variants> <gene_name>")
}

struct_variants <- args[1] # r"(C:\Users\qp241615\OneDrive - Queen Mary, University of London\Documents\4. Projects\1. DHX34\results\DDX41\DDX41_structural_variants.rds)"
gene_name       <- args[2] # "DDX41"


# Message validation files
cat("Validating input files for gene:", gene_name, "\n")
missing_files <- c()
if (!file.exists(struct_variants)) missing_files <- c(missing_files, paste("Structural variants:", struct_variants))

if (length(missing_files) > 0) {
  cat("ERROR: Missing files detected:\n", paste(missing_files, collapse = "\n"), "\n")
}
cat("All input files validated successfully\n")



# Libraries  -------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(circlize)
library(paletteer)


# Settings ----------------------------------------------------------------------
set.seed(23)

# Palettete ---------------------------------------------------------------------
my_pal <- c("#df4e50ff", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
            "#FFFF33", "#A65628", "#df75adff", "#844db8ff", "#66C2A5")

# Load Data ---------------------------------------------------------------------
## Structural Variants
variants_table <- readRDS(struct_variants)


# Plots Functions ---------------------------------------------------------
generate_barplot <- function(data, gene, variant_orig) {
  pdf(paste0(gene, "_", variant_orig, "_Frec_SVtype.pdf"), width = 10, height = 8)
  print(
    ggplot(data %>% count(INFO_SVTYPE), aes(x = INFO_SVTYPE, y = n, fill = INFO_SVTYPE)) +
      geom_col() +
      geom_text(aes(label = n), vjust = -0.5) +
      theme_minimal() +
      labs(title = paste0("Frequency of ", variant_orig, " variants for ", gene),
           x = "Type of Structural Variant", y = "Number of Structural Variants") +
      scale_fill_manual(values = my_pal) +
      theme(legend.position = "none"))
  dev.off()
}

generate_circos_plot <- function(data, gene, variant_orig) {
  pdf(paste0(gene, "_", variant_orig, "_Circle_Translocation.pdf"), width = 8, height = 8)
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
  pdf(paste0(gene, "_", variant_orig, "_CNV_Size_CN.pdf"), width = 10, height = 8)
  print(
    ggplot(data, aes(x = log10(SV_size), y = CN_CNVonly, color = CN_status)) +
      geom_point(alpha = 0.7, size = 3) +
      theme_minimal() +
      labs(
        title = paste0("Size vs Copy Number of ", variant_orig, " CNVs of ", gene),
        x = "log10(Size in bp)", y = "Copy Number"
      ) +
      scale_color_brewer(palette = "Set2", na.value = "grey50")
  )
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
}

cat("Plots generated successfully for gene:", gene_name, "\n")


