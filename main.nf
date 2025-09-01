#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { QC_SMALL_VARIANTS  } from './subworkflows/local/qc_small_variants/qc_small_variants.nf'
include { QC_STRUCT_VARIANTS } from './subworkflows/local/qc_struct_variants/qc_struct_variants.nf'

workflow {
        // Input validation
    if (!file(params.gene_list).exists()) {
        error "Gene list file not found: ${params.gene_list}"
    }
    if (!file(params.smallvar_annot_dir).exists()) {
        error "Small variants annotation directory not found: ${params.smallvar_annot_dir}"
    }
    if (!file(params.structvar_annot_dir).exists()) {
        error "Structural variants annotation directory not found: ${params.structvar_annot_dir}"
    }

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
    gene_list_ch = Channel
        .fromPath(params.gene_list)
        .splitText()
        .map { it.trim() }
        .filter { it }

    // Small Variant annotation files
    small_var_ch = Channel
        .fromPath("${params.smallvar_annot_dir}/*_annotated_variants.tsv")
        .filter { it.name.startsWith(params.genome_assembly) }

    // Structural Variant annotation files
    struct_var_ch = Channel
        .fromPath("${params.structvar_annot_dir}/*.tsv")
        .filter { it.name.contains("${params.genome_assembly}") }
        .filter { it.name.startsWith("CNV") || it.name.startsWith("SV") }
        

    // Run Small Variants QC workflow
    QC_SMALL_VARIANTS(
        small_var_ch,
        gene_list_ch,
        Channel.fromPath("${params.prot_dir}/*.gff"),
        Channel.fromPath("${params.exon_dir}/*.tsv")
    )

    // Run Structural Variants QC workflow
    QC_STRUCT_VARIANTS(
        struct_var_ch,
        gene_list_ch
    )

    workflow.onError = {
        log.error "Pipeline execution stopped with the following error: ${workflow.errorMessage}"
    }   
}
