// qc_small_variant.nf

include { cleanFormatSmallVar } from '../../../modules/local/R/smallVariants/cleanFormatSmallVar.nf'
include { extractSmallVarPartID } from '../../../modules/local/R/smallVariants/extractSmallVarPartID.nf'
include { plotQCsmallVar } from '../../../modules/local/R/smallVariants/plotQCsmallVar.nf'


workflow QC_SMALL_VARIANTS {
    take:
    small_var_ch    // channel: small variant files
    gene_list_ch    // channel: list of genes
    prot_files_ch   // channel: protein annotation files
    exon_files_ch   // channel: exon annotation files

    main:
    // Create versions channel
    ch_versions = Channel.empty()

    // Prepare input for cleanFormatSmallVar
    clean_input_ch = small_var_ch
        .combine(gene_list_ch)
        .filter { file, gene -> file.name.contains(gene) }
        .map { file, gene -> tuple(file, gene) }
    
    // RUN cleanFormatSmallVar
    cleanFormatSmallVar(clean_input_ch)  
            // out.clean_tsv: path(tsv_file) 
            // out.clean_rds: tuple(path(rds_file), val(gene_name))

    // RUN extractSmallVarPartID
    extractSmallVarPartID(cleanFormatSmallVar.out.clean_rds, params.labkey_main)
            // out.partID: path(partID_file)
            // out.partMet: path(partMetadata_file)

    // Prepare input for protein and exon files
    plot_input_ch = cleanFormatSmallVar.out.clean_rds
        .combine(prot_files_ch)
        .filter { _clean_table, gene_name, gff_file -> gff_file.name.contains(gene_name) }
        .combine(exon_files_ch)
        .filter { _clean_table, gene_name, _gff_file, exon_file -> exon_file.name.contains(gene_name) }
        .map { clean_table, gene_name, gff_file, exon_file -> tuple(clean_table, gene_name, gff_file, exon_file) }
    
    // RUN plotQCsmallVar
    plotQCsmallVar(plot_input_ch, extractSmallVarPartID.out.partMet)
            // out.plots: path(plot_file)

    // Collect versions
    ch_versions = ch_versions.mix(cleanFormatSmallVar.out.versions)
    ch_versions = ch_versions.mix(extractSmallVarPartID.out.versions)
    ch_versions = ch_versions.mix(plotQCsmallVar.out.versions)

    emit:
    clean_tsv   = cleanFormatSmallVar.out.clean_tsv     // File TSV
    clean_rds   = cleanFormatSmallVar.out.clean_rds     // tuple(path(rds_file), val(gene_name))
    partID      = extractSmallVarPartID.out.partID      // participant IDs .txt
    partMet     = extractSmallVarPartID.out.partMet     // participant metadata .tsv
    plots       = plotQCsmallVar.out[0]                 // PDFs
    versions    = ch_versions
}
