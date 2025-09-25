// qc_struct_variant.nf

include { cleanFormatStructVar_SV } from '../../../modules/local/R/structVariants/cleanFormatStructVar_SV.nf'
include { cleanFormatStructVar_CNV } from '../../../modules/local/R/structVariants/cleanFormatStructVar_CNV.nf'
include { cleanFormatStructVar } from '../../../modules/local/R/structVariants/cleanFormatStructVar.nf'
include { extractStructVarPartID } from '../../../modules/local/R/structVariants/extractStructVarPartID.nf'
include { plotQCstructVar } from '../../../modules/local/R/structVariants/plotQCstructVar.nf'
include { filterStructVar } from '../../../modules/local/R/structVariants/filterStructVar.nf'
include { plotFilteredStructVar } from '../../../modules/local/R/structVariants/plotFilteredStructVar.nf'

workflow QC_STRUCT_VARIANTS {
    take:
    struct_var_ch   // channel: structural variant files
    gene_list_ch    // channel: list of genes

    main:
    // Create versions channel
    ch_versions = Channel.empty()

    // Load reference annotation files
    prot_files_ch = Channel.fromPath("${params.prot_dir}/*.gff")
    exon_files_ch = Channel.fromPath("${params.exon_dir}/*.tsv")

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

    // Filter and prepare input channels for CNV and SV
    cnv_germline_ch = struct_var_ch
        .filter { it.name.startsWith("CNV") && it.name.contains("germline") }
    
    cnv_somatic_ch = struct_var_ch
        .filter { it.name.startsWith("CNV") && it.name.contains("somatic") }
    
    sv_germline_ch = struct_var_ch
        .filter { it.name.startsWith("SV") && it.name.contains("germline") }
    
    sv_somatic_ch = struct_var_ch
        .filter { it.name.startsWith("SV") && it.name.contains("somatic") }

    // Combine channels (germline and somatic) for CNV
    cnv_input_ch = cnv_germline_ch
        .combine(cnv_somatic_ch)
        .combine(gene_list_ch)
        .filter { germline_file, somatic_file, gene -> 
            germline_file.name.contains(gene) && somatic_file.name.contains(gene) 
        }
        .map { germline_file, somatic_file, gene -> 
            tuple(germline_file, somatic_file, gene) 
        }

    // Combine channels (germline and somatic) for SV
    sv_input_ch = sv_germline_ch
        .combine(sv_somatic_ch)
        .combine(gene_list_ch)
        .filter { germline_file, somatic_file, gene -> 
            germline_file.name.contains(gene) && somatic_file.name.contains(gene) 
        }
        .map { germline_file, somatic_file, gene -> 
            tuple(germline_file, somatic_file, gene) 
        }

    // RUN clean format for CNV and SV
    cleanFormatStructVar_CNV(cnv_input_ch)  
    cleanFormatStructVar_SV(sv_input_ch) 

    // Combine CNV and SV channels 
    cnv_clean_ch = cleanFormatStructVar_CNV.out.clean_cnv
        .map { file, gene -> tuple(gene, file) }  // Put gene first

    sv_clean_ch = cleanFormatStructVar_SV.out.clean_sv
        .map { file, gene -> tuple(gene, file) }   // Put gene first

    combined_struct_ch = cnv_clean_ch
        .join(sv_clean_ch, by: 0)  // Join by gene_name (now index 0)
        .map { gene_name, cnv_file, sv_file -> 
            tuple(cnv_file, sv_file, gene_name) 
        }

    // RUN cleanFormatStructVar
    cleanFormatStructVar(combined_struct_ch)

    // RUN extractStructVarPartID
    extractStructVarPartID(cleanFormatStructVar.out.clean_rds, params.labkey_main)
            // out.partID_txt: path(partID_file)
            // out.partMet_tsv: path(partMetadata_file)
            // out.partMet_rds: tuple(val(gene_name), path(partMetadata_rds_file))

    // Prepare input for filterStructVar - combine clean_rds with partMet
    cleanVar_Metadata_input_ch = cleanFormatStructVar.out.clean_rds
        .join(extractStructVarPartID.out.partMet_rds)  
        .map { gene_name, clean_table, part_met -> tuple(gene_name, clean_table, part_met) }

    // RUN plotQCstructVar
    plotQCstructVar(cleanVar_Metadata_input_ch)

    // RUN filterStructVar
    filterStructVar(cleanVar_Metadata_input_ch)
            // out.stats_csv: path(filtered_stats_file)
            // out.filtered_clean_tsv: path(filtered_tsv_file)
            // out.filtered_clean_rds: tuple(val(gene_name), path(filtered_rds_file))
    
    // Preparar input para plotFilteredStructVar 
    plot_filter_input_ch = filterStructVar.out.filtered_clean_rds
        .join(prot_files_with_gene_ch)
        .join(exon_files_with_gene_ch)
        .join(extractStructVarPartID.out.partMet_rds)
        .map { gene_name, filtered_table, prot_file, exon_file, part_met -> 
            tuple(gene_name, filtered_table, prot_file, exon_file, part_met) 
        }

    // RUN plotFilteredStructVar
    plotFilteredStructVar(plot_filter_input_ch)
            // out.filtered_plots: path(filtered_plot_file)

    // Collect versions
    ch_versions = ch_versions.mix(cleanFormatStructVar_CNV.out.versions)
    ch_versions = ch_versions.mix(cleanFormatStructVar_SV.out.versions)
    ch_versions = ch_versions.mix(cleanFormatStructVar.out.versions)
    ch_versions = ch_versions.mix(plotQCstructVar.out.versions) 
    if (params.enable_sql_queries) {
        ch_versions = ch_versions.mix(extractStructVarPartID.out.versions)
    }

    emit:
    clean_cnv_tables      = cnv_clean_ch.map { gene, file -> tuple(file, gene) } 
    clean_sv_tables       = sv_clean_ch.map { gene, file -> tuple(file, gene) }  
    partID_txt            = params.enable_sql_queries ? extractStructVarPartID.out.partID_txt : Channel.empty()
    partMet_tsv           = params.enable_sql_queries ? extractStructVarPartID.out.partMet_tsv : Channel.empty()
    partMet_rds           = params.enable_sql_queries ? extractStructVarPartID.out.partMet_rds : Channel.empty()
    filtered_clean_tsv    = filterStructVar.out.filtered_clean_tsv     
    plot_files            = plotQCstructVar.out.plots
    versions              = ch_versions
}



