// qc_small_variant.nf

include { cleanFormatSmallVar } from '../../../modules/local/R/smallVariants/cleanFormatSmallVar.nf'
include { extractSmallVarPartID } from '../../../modules/local/R/smallVariants/extractSmallVarPartID.nf'
include { plotQCsmallVar } from '../../../modules/local/R/smallVariants/plotQCsmallVar.nf'
include { filterSmallVar } from '../../../modules/local/R/smallVariants/filterSmallVar.nf'
include { plotFilteredSmallVar } from '../../../modules/local/R/smallVariants/plotFilteredSmallVar.nf'


workflow QC_SMALL_VARIANTS {
    take:
    small_var_ch    // channel: small variant files
    gene_list_ch    // channel: list of genes

    main:
    // Create versions channel
    ch_versions = Channel.empty()

    // Load reference annotation files
    prot_files_ch = Channel.fromPath("${params.prot_dir}/*.gff")
    exon_files_ch = Channel.fromPath("${params.exon_dir}/*.tsv")

    // Prepare input for cleanFormatSmallVar
    clean_input_ch = small_var_ch
        .combine(gene_list_ch)
        .filter { file, gene -> file.name.contains(gene) }
        .map { file, gene -> tuple(file, gene) }

    // RUN cleanFormatSmallVar
    cleanFormatSmallVar(clean_input_ch)  
            // out.clean_tsv: path(tsv_file) 
            // out.clean_rds: tuple(val(gene_name), path(rds_file))

    // RUN extractSmallVarPartID
    extractSmallVarPartID(cleanFormatSmallVar.out.clean_rds, params.labkey_main)
            // out.partID: path(partID_file)
            // out.partMet: path(partMetadata_file)
    
    // Prepare protein and exon files with gene_name as the key
    prot_files_with_gene_ch = prot_files_ch
        .map { file -> 
            def gene_name = file.baseName.replaceAll(/\.gff$/, '')
            tuple(gene_name, file)
        }

    exon_files_with_gene_ch = exon_files_ch
        .map { file -> 
            def gene_name = file.baseName.replaceAll(/\.tsv$/, '')
            tuple(gene_name, file)
        }

    // Prepare input for plotQCsmallVar
    plot_input_ch = cleanFormatSmallVar.out.clean_rds
        .join(prot_files_with_gene_ch)
        .join(exon_files_with_gene_ch)
        .combine(extractSmallVarPartID.out.partMet)
        .map { gene_name, clean_table, prot_file, exon_file, part_met -> 
            tuple(gene_name, clean_table, prot_file, exon_file, part_met) 
        }
    
    // RUN plotQCsmallVar
    plotQCsmallVar(plot_input_ch)
            // out.plots: path(plot_file)
    
    // Prepare input for filterSmallVar - combine clean_rds with partMet
    filter_input_ch = cleanFormatSmallVar.out.clean_rds
        .combine(extractSmallVarPartID.out.partMet)
        .map { gene_name, clean_table,  part_met -> tuple(gene_name, clean_table, part_met) }
    
    // RUN filterSmallVar
    filterSmallVar(filter_input_ch)
            // out.stats_csv: path(filtered_stats_file)
            // out.filtered_clean_tsv: path(filtered_tsv_file)
            // out.filtered_clean_rds: tuple(val(gene_name), path(filtered_rds_file))

    // Preparar input para plotFilteredSmallVar 
    plot_filter_input_ch = filterSmallVar.out.filtered_clean_rds
        .join(prot_files_with_gene_ch)
        .join(exon_files_with_gene_ch)
        .combine(extractSmallVarPartID.out.partMet)
        .map { gene_name, filtered_table, prot_file, exon_file, part_met -> 
            tuple(gene_name, filtered_table, prot_file, exon_file, part_met) 
        }

    // RUN plotFilteredSmallVar
    plotFilteredSmallVar(plot_filter_input_ch)
            // out.filtered_plots: path(filtered_plot_file)

    // Collect versions
    ch_versions = ch_versions.mix(cleanFormatSmallVar.out.versions)
    ch_versions = ch_versions.mix(extractSmallVarPartID.out.versions)
    ch_versions = ch_versions.mix(plotQCsmallVar.out.versions)
    ch_versions = ch_versions.mix(filterSmallVar.out.versions)
    ch_versions = ch_versions.mix(plotFilteredSmallVar.out.versions)

    emit:
    clean_tsv           = cleanFormatSmallVar.out.clean_tsv         // File TSV
    clean_rds           = cleanFormatSmallVar.out.clean_rds         // tuple(val(gene_name), path(rds_file)) 
    filtered_clean_tsv  = filterSmallVar.out.filtered_clean_tsv     // File TSV
    filtered_clean_rds  = filterSmallVar.out.filtered_clean_rds     // tuple(val(gene_name), path(rds_file))
    partID              = extractSmallVarPartID.out.partID          // participant IDs .txt
    partMet             = extractSmallVarPartID.out.partMet         // participant metadata .tsv
    plots               = plotQCsmallVar.out.plots                  // PDFs
    filtered_plots      = plotFilteredSmallVar.out.filtered_plots   // PDFs
    versions            = ch_versions
}
