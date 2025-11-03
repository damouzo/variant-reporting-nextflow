// qc_small_variant.nf

include { cleanFormatSmallVar } from '../../../modules/local/R/smallVariants/cleanFormatSmallVar.nf'
include { extractSmallVarPartID } from '../../../modules/local/R/smallVariants/extractSmallVarPartID.nf'
include { plotQCsmallVar } from '../../../modules/local/R/smallVariants/plotQCsmallVar.nf'
include { filterSmallVar } from '../../../modules/local/R/smallVariants/filterSmallVar.nf'
include { plotFilteredSmallVar } from '../../../modules/local/R/smallVariants/plotFilteredSmallVar.nf'
include { generateReport } from '../../../modules/local/quarto/generateReport.nf'


workflow QC_SMALL_VARIANTS {
    take:
    small_var_ch    // channel: small variant files
    gene_list_ch    // channel: list of genes

    main:
    // Create versions channel
    ch_versions = channel.empty()

    // Load reference annotation files
    prot_files_ch = channel.fromPath("${params.prot_dir}/*.gff")
    exon_files_ch = channel.fromPath("${params.exon_dir}/*.tsv")

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
        // out.partMet_rds: tuple(val(gene_name), path(partMet_rds_file))
        // out.partMet_tsv: path(partMet_tsv_file)
        // out.partID_txt: path(partID_txt_file)
    
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
        .join(extractSmallVarPartID.out.partMet_rds)
        .map { gene_name, clean_table, prot_file, exon_file, part_met -> 
            tuple(gene_name, clean_table, prot_file, exon_file, part_met) 
        }
    
    // RUN plotQCsmallVar
    plotQCsmallVar(plot_input_ch)
            // out.plots: path(plot_file)
    
    // Create unified filter channels that include both basic and custom filters
    unified_filter_input_ch = cleanFormatSmallVar.out.clean_rds
        .join(extractSmallVarPartID.out.partMet_rds)
        .flatMap { gene_name, clean_table, part_met ->
            def all_filters = []
            // Always add basic filtering (no custom filters)
            all_filters.add(tuple(gene_name, clean_table, part_met, 'filter_basic', [:]))
            // Add custom filters if they exist for this gene
            def gene_filters = params.custom_filters[gene_name]
            if (gene_filters) {
                gene_filters.each { filter_type, filter_config ->
                    if (filter_config.enabled) {
                        all_filters.add(tuple(gene_name, clean_table, part_met, filter_type, filter_config.filters))
                    }
                }
            }
            return all_filters
        }
    
    // RUN filterSmallVar
    filterSmallVar(unified_filter_input_ch)
            // out.stats_csv: path(filtered_stats_file)
            // out.filtered_clean_tsv: path(filtered_tsv_file)
            // out.filtered_clean_rds: tuple(val(gene_name), val(filter_type), path(filtered_rds_file))

    // Preparar input para plotFilteredSmallVar 
    // Combine filtered data with reference files and filtered metadata
    plot_filter_input_ch = filterSmallVar.out.filtered_clean_rds
        .join(filterSmallVar.out.filtered_metadata_rds, by: [0, 1])  // Join by gene_name and filter_type
        .combine(prot_files_with_gene_ch)
        .filter { gene_filter, _filter_type, _filtered_rds, _filtered_metadata, gene_prot, _prot_file -> gene_filter == gene_prot }
        .combine(exon_files_with_gene_ch)
        .filter { gene_filter, _filter_type, _filtered_rds, _filtered_metadata, _gene_prot, _prot_file, gene_exon, _exon_file -> gene_filter == gene_exon }
        .combine(extractSmallVarPartID.out.partMet_rds)
        .filter { gene_filter, _filter_type, _filtered_rds, _filtered_metadata, _gene_prot, _prot_file, _gene_exon, _exon_file, gene_meta_orig, _part_met -> gene_filter == gene_meta_orig }
        .map { gene_name, filter_type, filtered_table, filtered_metadata, _gene_prot, prot_file, _gene_exon, exon_file, _gene_meta_orig, part_met -> 
            tuple(gene_name, filter_type, filtered_table, prot_file, exon_file, part_met, filtered_metadata) 
        }

    // RUN plotFilteredSmallVar
    plotFilteredSmallVar(plot_filter_input_ch)
            // out.filtered_plots: path(filtered_plot_file)

    // Generate Small Variants QC Report
    template_file = channel.fromPath("${params.templates_dir}/quarto/small_variants_report.qmd")
    
    // Prepare QC plots keyed by gene
    qc_plots_keyed = plotQCsmallVar.out.plots
        .flatten()
        .map { file -> 
            def gene_name = file.baseName.split('_')[0]
            tuple(gene_name, file)
        }

    // Prepare filtered plots keyed by gene
    filtered_plots_keyed = plotFilteredSmallVar.out.filtered_plots
        .flatten()
        .map { file -> 
            def gene_name = file.baseName.split('_')[0]
            tuple(gene_name, file)
        }

    // Prepare filtered RDS files keyed by gene (group by gene, collect all filter types)
    filtered_rds_keyed = filterSmallVar.out.filtered_clean_rds
        .map { gene_name, _filter_type, rds_file -> tuple(gene_name, rds_file) }
        .groupTuple()

    // Prepare filtered statistics keyed by gene
    filtered_stats_keyed = filterSmallVar.out.stats_csv
        .flatten()
        .map { file -> 
            def gene_name = file.baseName.split('_')[0]
            tuple(gene_name, file)
        }

    // Prepare filtered TSV files keyed by gene
    filtered_tsv_keyed = filterSmallVar.out.filtered_clean_tsv
        .flatten()
        .map { file -> 
            def gene_name = file.baseName.split('_')[0]
            tuple(gene_name, file)
        }

    // Prepare participant files keyed by gene
    participant_txt_keyed = extractSmallVarPartID.out.partID_txt
        .flatten()
        .map { file -> 
            def gene_name = file.baseName.split('_')[0]
            tuple(gene_name, file)
        }

    participant_tsv_keyed = extractSmallVarPartID.out.partMet_tsv
        .flatten()
        .map { file -> 
            def gene_name = file.baseName.split('_')[0]
            tuple(gene_name, file)
        }

    // Prepare input files for the report - collect all files by gene
    report_files_ch = cleanFormatSmallVar.out.clean_rds
        .join(extractSmallVarPartID.out.partMet_rds)
        .join(
            qc_plots_keyed
                .groupTuple()
                .map { gene_name, plot_files -> tuple(gene_name, plot_files) }
        )
        .join(
            filtered_plots_keyed
                .groupTuple()
                .map { gene_name, filtered_plot_files -> tuple(gene_name, filtered_plot_files) }
        )
        .join(
            filtered_stats_keyed
                .groupTuple()
                .map { gene_name, stats_files -> tuple(gene_name, stats_files) }
        )
        .join(filtered_rds_keyed)
        .join(
            filtered_tsv_keyed
                .groupTuple()
                .map { gene_name, tsv_files -> tuple(gene_name, tsv_files) }
        )
        .join(
            participant_txt_keyed
                .groupTuple()
                .map { gene_name, txt_files -> tuple(gene_name, txt_files) }
        )
        .join(
            participant_tsv_keyed
                .groupTuple()
                .map { gene_name, part_tsv_files -> tuple(gene_name, part_tsv_files) }
        )
        .map { gene_name, clean_rds, part_met_rds, qc_plots, filtered_plots, stats_files, filtered_rds_files, filtered_tsv_files, part_txt_files, part_tsv_files ->
            // Combine all input files for this gene
            def input_files = [clean_rds, part_met_rds] + qc_plots + filtered_plots + stats_files + filtered_rds_files + filtered_tsv_files + part_txt_files + part_tsv_files
            def report_type = "small_variants_report"
            tuple(gene_name, report_type, input_files)
        }
        .combine(template_file)
        .map { gene_name, report_type, input_files, template ->
            tuple(gene_name, report_type, template, input_files) 
        }
        
    // Generate the Quarto report
    generateReport(report_files_ch)

    // Collect versions
    ch_versions = ch_versions.mix(cleanFormatSmallVar.out.versions)
    ch_versions = ch_versions.mix(extractSmallVarPartID.out.versions)
    ch_versions = ch_versions.mix(plotQCsmallVar.out.versions)
    ch_versions = ch_versions.mix(filterSmallVar.out.versions)
    ch_versions = ch_versions.mix(plotFilteredSmallVar.out.versions)
    ch_versions = ch_versions.mix(generateReport.out.versions)

    emit:
    clean_tsv               = cleanFormatSmallVar.out.clean_tsv         // File TSV
    clean_rds               = cleanFormatSmallVar.out.clean_rds         // tuple(val(gene_name), path(rds_file)) 
    filtered_clean_tsv      = filterSmallVar.out.filtered_clean_tsv     // File TSV (all filter types)
    filtered_clean_rds      = filterSmallVar.out.filtered_clean_rds     // tuple(val(gene_name), val(filter_type), path(rds_file))
    filtered_metadata_rds   = filterSmallVar.out.filtered_metadata_rds  // tuple(val(gene_name), val(filter_type), path(filtered_metadata_rds)) 
    filtered_stats_csv      = filterSmallVar.out.stats_csv              // Filter statistics (all filter types)
    partID_txt              = extractSmallVarPartID.out.partID_txt      // participant IDs .txt
    partMet_tsv             = extractSmallVarPartID.out.partMet_tsv     // participant metadata .tsv
    partMet_rds             = extractSmallVarPartID.out.partMet_rds     // participant metadata .rds
    plots                   = plotQCsmallVar.out.plots                  // PDFs
    filtered_plots          = plotFilteredSmallVar.out.filtered_plots   // PDFs (all filter types)
    html_reports            = generateReport.out.html_report            // HTML reports
    versions                = ch_versions
}
