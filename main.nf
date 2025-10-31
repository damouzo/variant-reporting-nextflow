#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { QC_SMALL_VARIANTS  } from './subworkflows/local/qc_small_variants/qc_small_variants.nf'
include { QC_STRUCT_VARIANTS } from './subworkflows/local/qc_struct_variants/qc_struct_variants.nf'
include { QC_BCI_PATIENTS_VARIANTS } from './subworkflows/local/qc_bci_patients_variants/qc_bci_patients_variants.nf'


workflow {
        // Input validation
    if (!file(params.gene_list).exists()) {error "Gene list not found: ${params.gene_list}"}
    if (!file(params.smallvar_annot_dir).exists()) {error "Small variants not found: ${params.smallvar_annot_dir}"}
    if (!file(params.structvar_annot_dir).exists()) {error "Structural variants not found: ${params.structvar_annot_dir}"}

    log.info """
    ===========================================
    ${params.project_code} : ${params.project_name}
    ===========================================
    Environment     : ${params.environment_type}
    ===========================================
    Gene List File  : ${params.gene_list}
    Results Dir     : ${params.results_dir}
    Work Dir        : ${workDir}
    ===========================================
    """

    // Gene list
    gene_list_ch = channel
        .fromPath(params.gene_list)
        .splitText()
        .map { line -> line.trim() }
        .filter { gene -> gene }

    // Small Variant annotation files
    small_var_ch = channel
        .fromPath("${params.smallvar_annot_dir}/*_annotated_variants.tsv")
        .filter { file -> file.name.startsWith(params.genome_assembly) }

    // Structural Variant annotation files
    struct_var_ch = channel
        .fromPath("${params.structvar_annot_dir}/*.tsv")
        .filter { file -> file.name.contains("${params.genome_assembly}") }
        .filter { file -> file.name.startsWith("CNV") || file.name.startsWith("SV") }

    // BCI patients variant files
    bci_patients_var_ch = channel
        .fromPath("${params.bci_patients_var_dir}/*.tsv")

    // Run Small Variants QC workflow
    QC_SMALL_VARIANTS(
        small_var_ch,
        gene_list_ch
    )

    // RUN Structural Variants QC workflow
    QC_STRUCT_VARIANTS(
        struct_var_ch,
        gene_list_ch
    )

    // RUN BCI samples QC workflow
    QC_BCI_PATIENTS_VARIANTS(
        bci_patients_var_ch,
        gene_list_ch,
        QC_SMALL_VARIANTS.out.clean_rds,
        QC_SMALL_VARIANTS.out.partMet_rds
    )
    
    workflow.onError = {
        log.error "Pipeline execution stopped with the following error: ${workflow.errorMessage}"
    }   
}
