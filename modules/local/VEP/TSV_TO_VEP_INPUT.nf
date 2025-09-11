// TSV_TO_VEP_INPUT.nf

process TSV_TO_VEP_INPUT {
    tag "${gene_name}"
    label 'r_process'

    input:
    tuple path(bci_patients_var_file), val(gene_name)

    output:
    tuple val(gene_name), path("${gene_name}_vep_input.txt"), emit: vep_input
    tuple val(gene_name), path("${gene_name}_original_data.tsv"), emit: original_data
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env Rscript
    
    library(readr)
    library(dplyr)
    
    # Read original data
    gene_name <- paste0("${gene_name}")
    df <- read_tsv("${bci_patients_var_file}")

    # Check for duplicates in cDNA_Change_HGVS_c.
    if (any(duplicated(df[, "cDNA_Change_HGVS_c."]))) {
        stop("Error: Duplicated values in 'cDNA_Change_HGVS_c.'. Please ensure unique identifiers.")
    }
    
    # Save original data for later merging
    write_tsv(df, "${gene_name}_original_data.tsv")
    
    # Prepare data for VEP
    vep_input_data <- df %>%
        filter(!is.na(cDNA_Change_HGVS_c.) & !is.na(Transcript)) %>%
        filter(cDNA_Change_HGVS_c. != "NA" & Transcript != "NA") %>%
        mutate(
            variant_id = paste(Patient_REF, row_number(), sep="_"), # Create unique identifier
            transcript_clean = sub("\\\\..*", "", Transcript), # Clean transcript (remove version if exists)
            vep_input_line = paste0(transcript_clean, ":", cDNA_Change_HGVS_c.) # Format for VEP: transcript:hgvs_c
        ) %>%
        select(variant_id, vep_input_line)
    
    # Write VEP input file (without headers)
    write_tsv(vep_input_data[,2], "${gene_name}_vep_input.txt", col_names = FALSE)

    cat("Prepared", nrow(vep_input_data), "variants of ", gene_name, " for VEP annotation\n")

    # Create versions file using R
    versions_content <- paste0(
        '"${task.process}":\\n',
        '    r-base: "', R.version.string, '"\\n',
        '    readr: "', packageVersion("readr"), '"\\n',
        '    dplyr: "', packageVersion("dplyr"), '"'
    )
    
    writeLines(versions_content, "versions.yml")

    """
}