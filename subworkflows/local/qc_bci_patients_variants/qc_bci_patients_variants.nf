// qc_bci_patients_variants.nf

include { TSV_TO_VEP_INPUT } from '../../../modules/local/VEP/TSV_TO_VEP_INPUT.nf'
include { VEP_ANNOTATE_TABULAR } from '../../../modules/local/VEP/VEP_ANNOTATE_TABULAR.nf'
include { cleanFormatBciPatientsVar } from '../../../modules/local/R/bciPatientsVariants/cleanFormatBciPatientsVar.nf'
include { plotQCBciPatientsVar } from '../../../modules/local/R/bciPatientsVariants/plotQCBciPatientsVar.nf'


workflow QC_BCI_PATIENTS_VARIANTS {
    take:
    bci_patients_var_ch     // channel: BCI patients variant files
    gene_list_ch            // channel: list of genes
    prot_files_ch           // channel: protein annotation files
    exon_files_ch           // channel: exon annotation files
    vep_cache               // path: VEP cache directory 

    main:
    // Create versions channel
    ch_versions = Channel.empty()

    // Prepare input for conversion to VEP format
    convert_input_ch = bci_patients_var_ch
        .combine(gene_list_ch)
        .filter { file, gene -> file.name.contains(gene) }
        .map { file, gene -> tuple(file, gene) }

    // Convert tabular data to VEP input format
    TSV_TO_VEP_INPUT(convert_input_ch)

    // Run VEP annotation
    VEP_ANNOTATE_TABULAR(
        TSV_TO_VEP_INPUT.out.vep_input,
        vep_cache,
    )

    // Merge VEP annotations with original data
    merge_input_ch = VEP_ANNOTATE_TABULAR.out.vep_output
        .join(TSV_TO_VEP_INPUT.out.original_data)

    // RUN cleanFormatBciPatients
    cleanFormatBciPatientsVar(merge_input_ch)
            // out.clean_tsv: path(tsv_file) 
            // out.clean_rds: path(rds_file)

    // Prepare input for protein and exon files
    plot_input_ch = cleanFormatBciPatientsVar.out.clean_rds
        .combine(prot_files_ch)
        .filter { _clean_table, gene_name, gff_file -> gff_file.name.contains(gene_name) }
        .combine(exon_files_ch)
        .filter { _clean_table, gene_name, _gff_file, exon_file -> exon_file.name.contains(gene_name) }
        .map { clean_table, gene_name, gff_file, exon_file -> tuple(clean_table, gene_name, gff_file, exon_file) }
    
    // RUN plotQCBciPatients
    plotQCBciPatientsVar(plot_input_ch)
            // out.plots: path(plot_file)

    // Collect versions
    ch_versions = ch_versions.mix(TSV_TO_VEP_INPUT.out.versions)
    ch_versions = ch_versions.mix(VEP_ANNOTATE_TABULAR.out.versions)
    ch_versions = ch_versions.mix(cleanFormatBciPatientsVar.out.versions)
    ch_versions = ch_versions.mix(plotQCBciPatientsVar.out.versions)

    emit:
    clean_tsv   = cleanFormatBciPatientsVar.out.clean_tsv     // File TSV
    clean_rds   = cleanFormatBciPatientsVar.out.clean_rds     // path(rds_file)
    plots       = plotQCBciPatientsVar.out.plots              // PDFs
    versions    = ch_versions
}
