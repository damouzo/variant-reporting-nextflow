#!/usr/bin/env Rscript
# Script: queryDiseaseRelated
# Queries 100kGP phenotype (HPO), ecds (SNOMED) and hes_apc (ICD-10) for disease terms
# and computes prevalence + enrichment against the whole 100kGP programme.

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
partID_txt  <- args[1]   # filtered participants .txt
gene_name   <- args[2]   # "DDX41"
filter_type <- args[3]   # "filter_basic"
labkey_main <- args[4]   # "/main-programme/main-programme_v19_2024-10-31"
terms_arg   <- args[5]   # "Anemia,271737000,HP:0001903"

# Libraries --------------------------------------------------------------------
library(tidyverse)
library(Rlabkey)

# LabKey helper ----------------------------------------------------------------
labkey.setWafEncoding(FALSE) # config for R>3
labkey_to_df <- function(sql_query, database, maxrows = 1000000) {
  labkey.setDefaults(baseUrl = "https://labkey-embassy.gel.zone/labkey/")
  labkey.executeSql(folderPath = database, schemaName = "lists",
                    colNameOpt = "rname", sql = sql_query, maxRows = maxrows) %>%
    mutate(across(everything(), as.character))
}

# Classify term -> ontology ----------------------------------------------------
classify <- function(t) {
  if (grepl("^HP:\\d{7}$", t))  return("HPO")
  if (grepl("^\\d{6,}$", t))    return("SNOMED")
  if (grepl("^[A-Z]\\d{2}", t)) return("ICD10")
  return("TEXT")
}

terms <- trimws(unlist(strsplit(terms_arg, ",")))
terms <- terms[terms != ""]
term_df <- tibble(term = terms, ontology = vapply(terms, classify, character(1)))

# Filtered participants + 100kGP total -----------------------------------------
filtered <- unique(readLines(partID_txt))
filtered <- filtered[filtered != ""]

prog_total <- labkey_to_df(
  "SELECT COUNT(DISTINCT participant_id) AS n FROM participant_summary", labkey_main)$n %>%
  as.numeric()

# One query per term -----------------------------------------------------------
query_term <- function(term, ontology) {
  sql <- switch(ontology,
    HPO    = paste0("SELECT DISTINCT \"Participant Id\" AS participant_id FROM rare_diseases_participant_phenotype ",
                    "WHERE lower(\"Hpo Present\") = 'yes' AND \"Hpo Id\" = '", term, "'"),
    SNOMED = paste0("SELECT DISTINCT participant_id FROM ecds WHERE diagnosis_code_all LIKE '%", term, "%'"),
    ICD10  = paste0("SELECT DISTINCT participant_id FROM hes_apc WHERE diag_all LIKE '%", term, "%'"),
    TEXT   = paste0("SELECT DISTINCT \"Participant Id\" AS participant_id FROM rare_diseases_participant_phenotype ",
                    "WHERE lower(\"Normalised Hpo Term\") LIKE '%", tolower(term), "%'")
  )
  tryCatch(labkey_to_df(sql, labkey_main), error = function(e) tibble(participant_id = character()))
}

cat("Querying", nrow(term_df), "disease term(s) in 100kGP...\n")
all_pos <- map_dfr(seq_len(nrow(term_df)), function(i) {
  df <- query_term(term_df$term[i], term_df$ontology[i])
  if (nrow(df) == 0) return(tibble(participant_id = character(), term = character()))
  transmute(df, participant_id, term = term_df$term[i])
})

# Stats + enrichment vs programme ----------------------------------------------
empty_stat <- tibble(term = character(), ontology = character(), n_filtered_total = numeric(),
                     n_filtered_with = numeric(), pct_filtered = numeric(),
                     n_programme_total = numeric(), n_programme_with = numeric(),
                     pct_programme = numeric(), odds_ratio = numeric(),
                     ci_low = numeric(), ci_high = numeric(), p_value = numeric())

stats <- if (nrow(term_df) == 0) empty_stat else {
  map_dfr(seq_len(nrow(term_df)), function(i) {
    t <- term_df$term[i]; ont <- term_df$ontology[i]
    n_prog_with <- sum(all_pos$term == t)
    n_filt_with <- sum(all_pos$term == t & all_pos$participant_id %in% filtered)
    a <- n_filt_with
    b <- length(filtered) - n_filt_with
    c <- n_prog_with - n_filt_with
    d <- (prog_total - length(filtered)) - c
    ft <- tryCatch(fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)),
                   error = function(e) NULL)
    tibble(term = t, ontology = ont,
           n_filtered_total = length(filtered), n_filtered_with = n_filt_with,
           pct_filtered = round(100 * n_filt_with / length(filtered), 2),
           n_programme_total = prog_total, n_programme_with = n_prog_with,
           pct_programme = round(100 * n_prog_with / prog_total, 2),
           odds_ratio = if (is.null(ft)) NA else unname(ft$estimate),
           ci_low  = if (is.null(ft)) NA else ft$conf.int[1],
           ci_high = if (is.null(ft)) NA else ft$conf.int[2],
           p_value = if (is.null(ft)) NA else ft$p.value)
  }) %>% mutate(fdr = p.adjust(p_value, method = "BH"))
}

# Outputs ----------------------------------------------------------------------
write.csv(stats, file = paste0(gene_name, "_", filter_type, "_disease_related_stats.csv"), row.names = FALSE)

all_pos %>%
  mutate(in_filtered = participant_id %in% filtered) %>%
  write.table(file = paste0(gene_name, "_", filter_type, "_disease_related_participants.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

# Plots (masking <5 only here) -------------------------------------------------
mask_lab <- function(n, pct) ifelse(n >= 1 & n <= 4, "1-4*", paste0(round(pct), "%"))
mask_val <- function(n, pct) ifelse(n >= 1 & n <= 4, 1, pct)

pdf(paste0(gene_name, "_", filter_type, "_disease_related_plots.pdf"), width = 11, height = 6)

bar_df <- bind_rows(
  tibble(term = stats$term, set = "Filtered",  pct = stats$pct_filtered,  n = stats$n_filtered_with),
  tibble(term = stats$term, set = "Programme", pct = stats$pct_programme, n = stats$n_programme_with)
) %>%
  mutate(label = ifelse(set == "Filtered", mask_lab(n, pct), paste0(round(pct), "%")),
         height = ifelse(set == "Filtered", mask_val(n, pct), pct))

print(ggplot(bar_df, aes(x = term, y = height, fill = set)) +
  geom_col(position = position_dodge(), alpha = 0.85) +
  geom_text(aes(label = label), position = position_dodge(width = 0.9), vjust = -0.4, size = 3) +
  labs(x = "Disease term", y = "% of participants",
       title = paste0(gene_name, " ", filter_type, ": disease prevalence (Filtered vs 100kGP)"),
       subtitle = "* masked counts (1-4)") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)))

forest_df <- stats %>% filter(!is.na(odds_ratio), is.finite(ci_low), is.finite(ci_high))
if (nrow(forest_df) > 0) {
  print(ggplot(forest_df, aes(x = odds_ratio, y = term)) +
    geom_point(size = 3, color = "#2c7fb8") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    scale_x_log10() +
    labs(x = "Odds ratio (log scale)", y = "",
         title = paste0(gene_name, " ", filter_type, ": enrichment vs 100kGP")) +
    theme_minimal())
}
dev.off()

cat("Done disease-related analysis for", gene_name, filter_type, "\n")
