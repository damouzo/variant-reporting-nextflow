// qc_small_variant.nf

include { cleanFormatSmallVar } from '../../../modules/local/R/smallVariants/cleanFormatSmallVar.nf'
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

    // Prepare input for protein and exon files
    plot_input_ch = cleanFormatSmallVar.out.clean_rds
        .combine(prot_files_ch)
        .filter { clean_table, gene_name, gff_file -> gff_file.name.contains(gene_name) }
        .combine(exon_files_ch)
        .filter { clean_table, gene_name, gff_file, exon_file -> exon_file.name.contains(gene_name) }
        .map { clean_table, gene_name, gff_file, exon_file -> tuple(clean_table, gene_name, gff_file, exon_file) }
    
    // RUN plotQCsmallVar
    plotQCsmallVar(plot_input_ch)

    // Collect versions
    ch_versions = ch_versions.mix(cleanFormatSmallVar.out.versions)
    ch_versions = ch_versions.mix(plotQCsmallVar.out.versions)

    emit:
    clean_tsv   = cleanFormatSmallVar.out.clean_tsv     // File TSV
    clean_rds   = cleanFormatSmallVar.out.clean_rds     // tuple(path(rds_file), val(gene_name))
    plots       = plotQCsmallVar.out[0]                 // PDFs
    versions    = ch_versions
}
