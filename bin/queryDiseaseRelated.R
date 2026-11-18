#!/usr/bin/env Rscript
# Script: queryDiseaseRelated
# Queries 100kGP phenotype (HPO), ecds (SNOMED) and hes_apc (ICD-10) for disease terms
# and computes prevalence + enrichment against the whole 100kGP programme.
# Also builds descriptive comorbidity barplots (all terms present among the
# filtered participants) per ontology.

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
partID_txt  <- args[1]   # filtered participants .txt
gene_name   <- args[2]   # "DDX41"
filter_type <- args[3]   # "filter_basic"
labkey_main <- args[4]   # "/main-programme/main-programme_v19_2024-10-31"
terms_arg   <- args[5]   # "Anemia,271737000,HP:0001903,D50"
exclude_arg <- if (length(args) >= 6) args[6] else "Z,U,V,W,X,Y"  # ICD-10 chapters dropped from comorbidity

# Libraries --------------------------------------------------------------------
library(tidyverse)
library(Rlabkey)

# LabKey helpers ---------------------------------------------------------------
labkey.setWafEncoding(FALSE) # config for R>3
labkey_to_df <- function(sql_query, database, maxrows = 1000000) {
  labkey.setDefaults(baseUrl = "https://labkey-embassy.gel.zone/labkey/")
  labkey.executeSql(folderPath = database, schemaName = "lists",
                    colNameOpt = "rname", sql = sql_query, maxRows = maxrows) %>%
    mutate(across(everything(), as.character))
}

# Raw read (one row) used only to discover the real column names of a table
labkey_raw <- function(sql_query, database, maxrows = 1) {
  labkey.setDefaults(baseUrl = "https://labkey-embassy.gel.zone/labkey/")
  labkey.executeSql(folderPath = database, schemaName = "lists",
                    colNameOpt = "fieldname", sql = sql_query, maxRows = maxrows)
}

discover_columns <- function(table, database) {
  tryCatch({
    cols <- names(labkey_raw(paste0("SELECT * FROM ", table, " LIMIT 1"), database, 1))
    cat("  columns of", table, ":", paste(cols, collapse = " | "), "\n")
    cols
  }, error = function(e) {
    cat("  ! could not read columns of", table, ":", conditionMessage(e), "\n")
    character(0)
  })
}

pick_col <- function(cols, patterns) {
  for (p in patterns) {
    hit <- cols[grepl(p, cols, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
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

# ICD-10 chapters to exclude from the descriptive comorbidity analysis (they are
# not diseases: Z = factors, U = special, V-Y = external causes).
icd10_exclude <- trimws(unlist(strsplit(exclude_arg, ",")))
icd10_exclude <- icd10_exclude[icd10_exclude != ""]

# Discover the real column names of the phenotype table so the query works
# regardless of the exact naming used by the LabKey dataset.
cat("Inspecting LabKey tables...\n")
ph_cols <- discover_columns("rare_diseases_participant_phenotype", labkey_main)
ph_pid  <- pick_col(ph_cols, c("participant.*id"))
ph_hpo  <- pick_col(ph_cols, c("hpo.*id", "hpo.*code"))
ph_pres <- pick_col(ph_cols, c("hpo.*present", "present"))
ph_term <- pick_col(ph_cols, c("normalis.*term", "normaliz.*term", "hpo.*term", "term"))
ph_pid  <- if (!is.na(ph_pid))  ph_pid  else "Participant Id"
ph_hpo  <- if (!is.na(ph_hpo))  ph_hpo  else "Hpo Id"
ph_pres <- if (!is.na(ph_pres)) ph_pres else "Hpo Present"
ph_term <- if (!is.na(ph_term)) ph_term else "Normalised Hpo Term"
cat("  using phenotype columns:", ph_pid, "/", ph_hpo, "/", ph_pres, "/", ph_term, "\n")

# Spelling-tolerant match for free-text phenotype terms (US "anemia" vs UK "anaemia").
text_term_condition <- function(term, col) {
  vars <- unique(tolower(trimws(c(term, gsub("e", "ae", term), gsub("ae", "e", term)))))
  conds <- paste0("lower(\"", col, "\") LIKE '%", vars, "%'")
  paste0("(", paste(conds, collapse = " OR "), ")")
}

# Filtered participants + 100kGP total -----------------------------------------
filtered <- unique(readLines(partID_txt))
filtered <- filtered[filtered != ""]

prog_total <- labkey_to_df(
  "SELECT COUNT(DISTINCT participant_id) AS n FROM participant_summary", labkey_main)$n %>%
  as.numeric()

# One query per term -----------------------------------------------------------
query_term <- function(term, ontology) {
  sql <- switch(ontology,
    HPO    = paste0("SELECT DISTINCT \"", ph_pid, "\" AS participant_id FROM rare_diseases_participant_phenotype ",
                    "WHERE lower(\"", ph_pres, "\") IN ('yes','y','true','1','present') AND \"", ph_hpo, "\" = '", term, "'"),
    SNOMED = paste0("SELECT DISTINCT participant_id FROM ecds WHERE diagnosis_code_all LIKE '%", term, "%'"),
    ICD10  = paste0("SELECT DISTINCT participant_id FROM hes_apc WHERE diag_all LIKE '%", term, "%'"),
    TEXT   = paste0("SELECT DISTINCT \"", ph_pid, "\" AS participant_id FROM rare_diseases_participant_phenotype ",
                    "WHERE ", text_term_condition(term, ph_term))
  )
  tryCatch(labkey_to_df(sql, labkey_main),
           error = function(e) { cat("    ! query error:", conditionMessage(e), "\n"); tibble(participant_id = character()) })
}

cat("Querying", nrow(term_df), "disease term(s) in 100kGP...\n")
all_pos <- map_dfr(seq_len(nrow(term_df)), function(i) {
  df <- query_term(term_df$term[i], term_df$ontology[i])
  cat("  -", term_df$term[i], "(", term_df$ontology[i], "):", nrow(df), "participant(s)\n")
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

# Helpers: masking + significance ----------------------------------------------
p_to_stars <- function(p) dplyr::case_when(
  is.na(p)  ~ "",
  p < 0.001 ~ "***",
  p < 0.01  ~ "**",
  p < 0.05  ~ "*",
  TRUE      ~ "ns"
)

# Masked cells (1-4 participants) keep the label but get a fixed, tiny height so
# they are not visually mistaken for a real prevalence.
mask_lab <- function(n, pct) ifelse(n >= 1 & n <= 4, "1-4*", paste0(sprintf("%.1f", pct), "%"))
mask_val <- function(n, pct) ifelse(n >= 1 & n <= 4, 0.1, pct)

# ---------------------------------------------------------------------------
# Plot 1: prevalence barplot (Filtered vs 100kGP) + Fisher significance
# ---------------------------------------------------------------------------
term_levels <- stats$term
K_terms <- length(term_levels)

bar_df <- bind_rows(
  tibble(term = stats$term, set = "Filtered",  pct = stats$pct_filtered,  n = stats$n_filtered_with),
  tibble(term = stats$term, set = "Programme", pct = stats$pct_programme, n = stats$n_programme_with)
) %>%
  mutate(i = match(term, term_levels),
         label = ifelse(set == "Filtered", mask_lab(n, pct), paste0(sprintf("%.1f", pct), "%")),
         height = ifelse(set == "Filtered", mask_val(n, pct), pct))

has_masked <- any(stats$n_filtered_with >= 1 & stats$n_filtered_with <= 4)
subtitle_txt <- if (has_masked) "* masked counts (1-4)" else NULL

p_bar <- ggplot(bar_df, aes(x = i, y = height, fill = set)) +
  geom_col(position = position_dodge(width = 0.9), alpha = 0.85) +
  geom_text(aes(label = label), position = position_dodge(width = 0.9), vjust = -0.4, size = 3) +
  scale_x_continuous(breaks = seq_len(K_terms), labels = term_levels) +
  labs(x = "Disease term", y = "% of participants",
       title = paste0(gene_name, " ", filter_type, ": disease prevalence (Filtered vs 100kGP)"),
       subtitle = subtitle_txt) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Significance brackets (Fisher vs rest of programme), annotated with FDR stars.
# Bars are dodged, so the bracket connects i - 0.225 and i + 0.225.
stars_df <- stats %>%
  mutate(stars = p_to_stars(fdr)) %>%
  filter(stars %in% c("*", "**", "***")) %>%
  mutate(i = match(term, term_levels))

if (nrow(stars_df) > 0) {
  ymax_df <- bar_df %>% group_by(i) %>% summarise(ymax = max(height), .groups = "drop")
  off <- max(0.15, 0.06 * max(bar_df$height))
  stars_df <- stars_df %>% left_join(ymax_df, by = "i") %>% mutate(y = ymax + off)
  p_bar <- p_bar +
    geom_segment(data = stars_df, aes(x = i - 0.225, xend = i + 0.225, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = stars_df, aes(x = i - 0.225, xend = i - 0.225, y = y, yend = y - off / 2),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = stars_df, aes(x = i + 0.225, xend = i + 0.225, y = y, yend = y - off / 2),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(data = stars_df, aes(x = i, y = y + off / 2, label = stars),
              inherit.aes = FALSE, size = 4)
}

# ---------------------------------------------------------------------------
# Plot 2: forest plot of enrichment
# ---------------------------------------------------------------------------
forest_df <- stats %>% filter(!is.na(odds_ratio), is.finite(ci_low), is.finite(ci_high))

pdf(paste0(gene_name, "_", filter_type, "_disease_related_plots.pdf"), width = 11, height = 6)
print(p_bar)
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

# ---------------------------------------------------------------------------
# Plot 3: comorbidity (descriptive, no statistics)
# All disease terms present among the filtered participants, per ontology.
# 100% = number of filtered participants. Terms with 1-4 participants are
# merged into "Other / small counts masked" (counted as unique participants).
# ---------------------------------------------------------------------------
query_all_terms <- function(ids, ontology) {
  id_list <- paste0("'", ids, "'", collapse = ", ")
  sql <- switch(ontology,
    SNOMED = paste0("SELECT DISTINCT participant_id, diagnosis_code_all AS codes FROM ecds ",
                    "WHERE participant_id IN (", id_list, ")"),
    ICD10  = paste0("SELECT DISTINCT participant_id, diag_all AS codes FROM hes_apc ",
                    "WHERE participant_id IN (", id_list, ")"),
    HPO    = paste0("SELECT DISTINCT \"", ph_pid, "\" AS participant_id, \"", ph_hpo, "\" AS codes ",
                    "FROM rare_diseases_participant_phenotype ",
                    "WHERE lower(\"", ph_pres, "\") IN ('yes','y','true','1','present') AND \"", ph_pid, "\" IN (", id_list, ")")
  )
  tryCatch(labkey_to_df(sql, labkey_main),
           error = function(e) { cat("    ! comorbidity query error:", conditionMessage(e), "\n")
                                 tibble(participant_id = character(), codes = character()) })
}

explode_codes <- function(df) {
  if (nrow(df) == 0) return(tibble(participant_id = character(), term = character()))
  df %>%
    separate_rows(codes, sep = "[|,;]") %>%
    transmute(participant_id, term = str_trim(codes)) %>%
    filter(!is.na(term), term != "", term != "NA") %>%
    distinct(participant_id, term)
}

make_comorbidity_plot <- function(ontology, nice_label) {
  long <- explode_codes(query_all_terms(filtered, ontology))
  cat("Comorbidity", ontology, ":", nrow(long), "participant-term pair(s)\n")
  if (nrow(long) == 0) return(invisible(NULL))

  if (ontology == "ICD10") {
    # group sub-codes to the 3-character ICD-10 category and drop non-disease
    # chapters (Z/U/V-Y)
    long <- long %>%
      mutate(term = str_sub(term, 1, 3)) %>%
      filter(!str_sub(term, 1, 1) %in% icd10_exclude)
  }

  n_total  <- length(filtered)
  per_term <- long %>% count(term, sort = TRUE)
  visible  <- per_term %>% filter(n >= 5) %>% head(20)
  if (nrow(visible) == 0) return(invisible(NULL))

  kept    <- visible$term
  n_kept  <- long %>% filter(term %in% kept) %>% distinct(participant_id) %>% nrow()
  other_n <- max(0, n_total - n_kept)

  plot_df <- bind_rows(
    visible %>% transmute(term, n, pct = 100 * n / n_total, grp = "term"),
    tibble(term = "Other / small counts masked", n = other_n,
           pct = 100 * other_n / n_total, grp = "other")
  ) %>%
    mutate(term = fct_reorder(term, pct))

  p <- ggplot(plot_df, aes(x = pct, y = term, fill = grp)) +
    geom_col(alpha = 0.85) +
    geom_text(aes(label = paste0(sprintf("%.1f", pct), "%")), hjust = -0.1, size = 3) +
    scale_fill_manual(values = c(term = "#2c7fb8", other = "#b8b8b8"), guide = "none") +
    labs(x = paste0("% of ", n_total, " filtered participants"), y = "",
         title = paste0(gene_name, " ", filter_type, ": ", nice_label, " terms (descriptive)"),
         subtitle = "Terms present among filtered participants; counts 1-4 merged") +
    expand_limits(x = max(plot_df$pct) * 1.15) +
    theme_minimal()

  ggsave(paste0(gene_name, "_", filter_type, "_disease_related_comorbidity_", ontology, ".pdf"),
         plot = p, width = 11, height = max(5, 0.35 * nrow(plot_df) + 2))
  invisible(NULL)
}

cat("Building comorbidity plots...\n")
make_comorbidity_plot("SNOMED", "SNOMED")
make_comorbidity_plot("ICD10", "ICD-10")
make_comorbidity_plot("HPO", "HPO")

cat("Done disease-related analysis for", gene_name, filter_type, "\n")
