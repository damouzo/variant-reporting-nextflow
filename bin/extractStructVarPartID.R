#!/usr/bin/env Rscript
# Script: extractStructVarPartID

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
annotation_file <- args[1]  # "C:/Users/qp241615/OneDrive - Queen Mary, University of London/Documents/4. Projects/1. DHX34/data/raw_data/WGS_Variants/smallVar/GRCh38_DDX41_annotated_variants.tsv"
gene_name <- args[2]        # "DDX41"
labkey_main <- args[3]      # "/main-programme/main-programme_v19_2024-10-31"


# Libraries  -------------------------------------------------------------------
library(tidyverse)
library(Rlabkey)

# Settings ----------------------------------------------------------------------
set.seed(23)

# Load gene annotated variant file ----------------------------------------------
variant_table <- readRDS(annotation_file)


# SQL to Labkey -----------------------------------------------------------------
    # Extract samples of variants 
sample_vector <- unique(variant_table$SAMPLE)
sample_vector <- sample_vector[!is.na(sample_vector) & sample_vector != ""] 
all_samples <- unique(sample_vector[sample_vector != ""])


    # Config Labkey
labkey.setWafEncoding(FALSE) #config for R>3
sample_list_sql <- paste0("'", all_samples, "'", collapse = ", ")

sql <- paste0("
  SELECT
    sr.participant_id,
    sr.plate_key,
    sr.lab_sample_id,
    sr.genome_build,
    sr.path,
    ps.*
  FROM
    sequencing_report sr
  JOIN
    participant_summary ps ON 
    sr.participant_id = ps.participant_id
  WHERE
    plate_key IN (", sample_list_sql, ")
")

    # SQL Query
labkey_to_df <- function(sql_query, database, maxrows){
  labkey.setDefaults(baseUrl = "https://labkey-embassy.gel.zone/labkey/")
  labkey.executeSql(folderPath = database,
                    schemaName = "lists",
                    colNameOpt = "rname",
                    sql = sql_query,
                    maxRows = maxrows) %>%
    mutate(across(everything(), as.character))
}

query <- labkey_to_df(sql, labkey_main, 1000000)

    # Condense by Participant
participant_metadata <- query %>%
    group_by(participant_id) %>%
    summarise(across(everything(), ~ paste(unique(.[!is.na(.)]), collapse = ", ")),
              .groups = "drop")



# Output ------------------------------------------------------------------------
cat("/nGlimpse of metadata and ID files of participants")

# Save participants metadata
write.table(participant_metadata, file = paste0(gene_name, "_structural_variants_participantMetadata.tsv"), 
      sep = "\t", row.names = FALSE, quote = FALSE)

# Save participants metadata as RDS
saveRDS(participant_metadata, file = paste0(gene_name, "_structural_variants_participantMetadata.rds"))

# Save all participant IDs
all_participants <- unique(query$participant_id)
writeLines(all_participants, paste0(gene_name, "_structural_variants_participantID.txt"))
