#!/usr/bin/env Rscript

# Load required packages, installing if necessary
required_packages <- c("ggplot2", "dplyr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# --- Get command-line arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("script.R <gene_name> <clean_table> <prot_file> <exon_file>")
}

clean_table_file <- args[1]
gene_name        <- args[2]
prot_file        <- args[3]
exon_file        <- args[4]

cat("QC plotting for gene:", gene_name, "\n")
cat("clean_table_file:", clean_table_file, "\n")
cat("prot_file:", prot_file, "\n")
cat("exon_file:", exon_file, "\n")

# --- Datos de ejemplo ultrabásicos ---
df_bar <- data.frame(
  categoria = c("A", "B", "C"),
  valor = c(3, 7, 2)
)

df_scatter <- data.frame(
  x = 1:10,
  y = c(2, 5, 3, 6, 8, 7, 5, 9, 4, 6)
)

# --- Crear carpeta de salida ---
out_dir <- file.path("plots", "small_variants")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# --- Plot 1: Barplot básico ---
p1 <- ggplot(df_bar, aes(x = categoria, y = valor, fill = categoria)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Barplot básico", x = "Categoría", y = "Valor")

pdf(
  file.path(out_dir, paste0(gene_name, "_barplot.pdf")),
  width = 5, height = 4
)
print(p1)
dev.off()

# --- Plot 2: Scatterplot básico ---
p2 <- ggplot(df_scatter, aes(x = x, y = y)) +
  geom_point(size = 2, color = "blue") +
  theme_minimal() +
  labs(title = "Scatterplot básico", x = "X", y = "Y")

pdf(
  file.path(out_dir, paste0(gene_name, "_scatterplot.pdf")),
  width = 5, height = 4
)
print(p2)
dev.off()

cat("Plots básicos guardados en:", out_dir, "\n")
