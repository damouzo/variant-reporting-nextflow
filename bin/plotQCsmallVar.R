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
library(ggrepel)
library(Rlabkey)
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
densities <- variants_table %>% 
  group_by(IMPACT_annotation) %>%
  group_map(~ density(.x$POS_variant)$y)

max_density <- max(unlist(densities))

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

    
## Lollipop -------
# Check if variants have canonical annotation
if("CANONICAL_annotation" %in% colnames(variants_table)) {
  for (canonical_var in unique(variants_table$CANONICAL_annotation)) {
    
    if (canonical_var == "YES") {canonical_label <- "Canonical"} else {canonical_label <- "Non-canonical"}
  
    df_lollipop <- variants_table %>% 
      dplyr::filter(!is.na(Protein_pos_start), !is.na(NS_variant), NS_variant>0, 
                    IMPACT_annotation %in% c("HIGH"),CANONICAL_annotation == canonical_var)
    
    top_variants <- df_lollipop %>% group_by(alignment) %>%
      slice_max(order_by = NS_variant, n = 10, with_ties = F) %>%  ungroup()
    
    # Obtain protein annotation files are available
    protein_domains <- as.data.frame(protein_info) %>%
      filter(type %in% c("Domain", "Motif", "Region", "Zinc finger")) %>%
      mutate(type = factor(type, levels = c("Domain", "Region", "Motif")),
        ymin = case_when(type == "Domain" ~ -5,
                          type == "Region" ~ -10,
                          type == "Motif"  ~ -15, TRUE ~ -10),
        ymax = ymin + 3,
        ytext = ymin - 2  
      )
    
    
    pdf(paste0(gene_name,"_", canonical_label,"_Lollipop.pdf"), width=14, height=8)
    p <- ggplot(df_lollipop, aes(x = Protein_pos_start))
    
    # Add domains if available
    if(nrow(protein_domains) > 0) {
      p <- p +
        geom_rect(data = protein_domains,aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = type),
                  inherit.aes = FALSE, color = "black", alpha = 0.4) +
        geom_text(data = protein_domains,aes(x = (start + end)/2, y = ytext, label = Note),
                  inherit.aes = FALSE, size = 3, angle = 0, hjust = 0.5)
    }
    
    # Add variants
    p <- p +
      geom_segment(aes(xend = Protein_pos_start, y = 0, yend = NS_variant),
                   color = "grey60", size = 0.7) +
      geom_point(aes(y = NS_variant, color = Consequence_annotation, shape = BIOTYPE_annotation),
                 size = 5, stroke = 1, fill = "white") +
      geom_text_repel(data = top_variants,aes(x=Protein_pos_start, y=NS_variant, label = paste0("Max AF ", MAX_AF_annotation)),
                      size = 3.5, min.segment.length = 0, max.overlaps = Inf, box.padding = 0.4) +
      scale_color_manual(values = my_pal, name = "Consequence") +
      scale_shape(name = "Biotype") +
      facet_wrap(~ alignment, ncol = 1, strip.position = "top") +
      labs(title = "Lollipop plot over protein",
           subtitle = paste0(canonical_label, " High IMPACT variants, with prot_pos and NS_variant"),
           x = "Protein position", y = "Number of samples") +
      theme_minimal(base_size = 14) +
      theme(panel.grid.major.y = element_line(color = "grey90"),panel.grid.major.x = element_blank(),
            axis.text.y = element_text(size = 10),strip.text = element_text(face = "bold", size = 14))
    
    print(p)
    dev.off() 
    
    
    ## Lollipop 2.0 ----------------------------------------------------------------
    # Only create Lollipop 2.0 if protein annotation is available
    
    # Domains
    feature2plot <- c("Domain", "Motif", "Region", "Zinc finger")
    features <- protein_info[protein_info$type %in% feature2plot, ]
          
    # Prepare colors and labels
    features$Note <- trimws(features$Note)
    note_colors <- setNames(my_pal[1:length(unique(features$Note))], unique(features$Note))
    features$fill <- note_colors[features$Note]
    names(features) <- features$Note
      
    
    # Variants
    sample.gr <- GRanges(seqnames = paste0(unique(seqnames(features))), 
                          ranges = IRanges(start = df_lollipop$Protein_pos_start, width = 1, 
                                          names = df_lollipop$LabelVarPlot))
    
    # Prepare colors and labels
    sample.gr$score <- df_lollipop$NS_variant
    
    sample.gr$Consequence_annotation <- df_lollipop$Consequence_annotation
    sample.gr$color <- setNames(pastel_colors, unique(sample.gr$Consequence_annotation))[sample.gr$Consequence_annotation]
    legends <- list(labels=unique(sample.gr$Consequence_annotation), fill=unique(sample.gr$color))
    
    sample.gr$shape <- ifelse(grepl("^rs", names(sample.gr)), "circle", "diamond")
    
    sample.gr.rot <- sample.gr
    sample.gr.rot$label.parameter.rot <- 45
    if (canonical_var != "YES") {names(sample.gr.rot) <- NULL}
  
    # Plot Per Se       
    pdf(paste0(gene_name,"_", canonical_label,"_Lollipop2.0.pdf"))
    lolliplot(sample.gr.rot, features, legend=legends,ylab="Num. of Samples")
    grid.text(paste0("High ",canonical_label," Variants of ", gene_name), x=.5, y=.98,gp=gpar(cex=1.5, fontface="bold"))
    dev.off() 
    
  }
} else {
  cat("No CANONICAL_annotation column found. Skipping lollipop plots.\n")
}


cat("QC plotting for gene:", gene_name, "\n")
cat("Plots generated successfully!\n")
