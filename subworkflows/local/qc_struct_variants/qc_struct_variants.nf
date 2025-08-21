// qc_struct_variant.nf

include { cleanFormatStructVar_SV } from '../../../modules/local/R/structVariants/cleanFormatStructVar_SV.nf'
include { cleanFormatStructVar_CNV } from '../../../modules/local/R/structVariants/cleanFormatStructVar_CNV.nf'
include { cleanFormatStructVar } from '../../../modules/local/R/structVariants/cleanFormatStructVar.nf'
include { plotQCstructVar } from '../../../modules/local/R/structVariants/plotQCstructVar.nf'


workflow QC_STRUCT_VARIANTS {
    take:
    struct_var_ch   // channel: structural variant files
    gene_list_ch    // channel: list of genes

    main:
    // Create versions channel
    ch_versions = Channel.empty()

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
    cnv_clean_ch = cleanFormatStructVar_CNV.out[0]
        .map { file, gene -> tuple(gene, file) }  // Put gene first
    
    sv_clean_ch = cleanFormatStructVar_SV.out[0]
        .map { file, gene -> tuple(gene, file) }   // Put gene first

    combined_struct_ch = cnv_clean_ch
        .join(sv_clean_ch, by: 0)  // Join by gene_name (now index 0)
        .map { gene_name, cnv_file, sv_file -> 
            tuple(cnv_file, sv_file, gene_name) 
        }

    // RUN cleanFormatStructVar
    cleanFormatStructVar(combined_struct_ch)

    // RUN plotQCstructVar
    plotQCstructVar(cleanFormatStructVar.out.struct_rds)

    // Collect versions
    ch_versions = ch_versions.mix(cleanFormatStructVar_CNV.out.versions)
    ch_versions = ch_versions.mix(cleanFormatStructVar_SV.out.versions)
    ch_versions = ch_versions.mix(cleanFormatStructVar.out.versions)

    emit:
    clean_cnv_tables = cnv_clean_ch.map { gene, file -> tuple(file, gene) }  // Return to original format
    clean_sv_tables  = sv_clean_ch.map { gene, file -> tuple(file, gene) }   // Return to original format
    final_tables     = cleanFormatStructVar.out[0]
    versions         = ch_versions
}

