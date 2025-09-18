// qc_bci_patients_variants.nf

include { TSV_TO_VEP_INPUT } from '../../../modules/local/VEP/TSV_TO_VEP_INPUT.nf'
include { VEP_ANNOTATE_TABULAR } from '../../../modules/local/VEP/VEP_ANNOTATE_TABULAR.nf'
include { cleanFormatBciPatientsVar } from '../../../modules/local/R/bciPatientsVariants/cleanFormatBciPatientsVar.nf'
include { compareBCIwithGEvar } from '../../../modules/local/R/bciPatientsVariants/compareBCIwithGEvar.nf'
include { plotQCBciPatientsVar } from '../../../modules/local/R/bciPatientsVariants/plotQCBciPatientsVar.nf'
include { plotMetadataBCIvsGE } from '../../../modules/local/R/bciPatientsVariants/plotMetadataBCIvsGE.nf'



workflow QC_BCI_PATIENTS_VARIANTS {
    take:
    bci_patients_var_ch     // channel: BCI patients variant files
    gene_list_ch            // channel: list of genes
    prot_files_ch           // channel: protein annotation files
    exon_files_ch           // channel: exon annotation files
    clean_variants_ch       // channel: GE small variants from SmallVar subworkflow
    part_metadata_ch        // channel: participant metadata from SmallVar subworkflow

    main:
    // Create versions channel
    ch_versions = Channel.empty()

    
    // Conditional logic based on environment capabilities
    if (params.enable_vep_annotation) {
        // VEP annotation path (for personal_pc and hpc_apocrita)
        log.info "Using VEP annotation workflow for BCI patients variants"
        
        // Prepare input for conversion to VEP format
        convert_input_ch = bci_patients_var_ch
            .combine(gene_list_ch)
            .filter { file, gene -> file.name.contains(gene) }
            .map { file, gene -> tuple(file, gene) }

        // Convert tabular data to VEP input format
        TSV_TO_VEP_INPUT(convert_input_ch)

        // RUN VEP annotation
        VEP_ANNOTATE_TABULAR(TSV_TO_VEP_INPUT.out.vep_input)

        // Merge VEP annotations with original data
        merge_input_ch = VEP_ANNOTATE_TABULAR.out.vep_output
            .join(TSV_TO_VEP_INPUT.out.original_data)

        // RUN Clean and format BCI patients variants
        cleanFormatBciPatientsVar(merge_input_ch)

        // Collect VEP versions
        ch_versions = ch_versions.mix(TSV_TO_VEP_INPUT.out.versions)
        ch_versions = ch_versions.mix(VEP_ANNOTATE_TABULAR.out.versions)

    } else {
        // Pre-annotated data path (for ge_pc)
        log.info "Using pre-annotated BCI patients variants from ${params.bci_annotated_dir}"
        
        // Use pre-annotated files
        annotated_files_ch = Channel
            .fromPath("${params.bci_annotated_dir}/*_vep_output.txt")
            .combine(gene_list_ch)
            .filter { file, gene -> file.name.contains(gene) }
            .map { file, gene -> tuple(gene, file) }

        // Create channels for original files  
        original_files_ch = Channel
            .fromPath("${params.bci_annotated_dir}/*_patients_var_data.tsv")
            .combine(gene_list_ch)
            .filter { file, gene -> file.name.contains(gene) }
            .map { file, gene -> tuple(gene, file) }

        // Combine annotated and original files by gene
        pre_annotated_ch = annotated_files_ch
            .join(original_files_ch)  // Une por gene: [gene, annotated_file, original_file]

        // Debug: ver qué archivos están llegando
        pre_annotated_ch.view { "DEBUG pre_annotated_ch: gene=${it[0]}, annotated=${it[1].name}, original=${it[2].name}" }      

        // Clean and format BCI patients variants (same process, different input)
        cleanFormatBciPatientsVar(pre_annotated_ch)
    }
    
    // Combine BCI clean data with clean variants from the first subworkflow
    comparison_input_ch = cleanFormatBciPatientsVar.out.clean_rds
        .combine(clean_variants_ch)
        .filter { _bci_rds, bci_gene, _ge_rds, ge_gene -> 
            bci_gene == ge_gene 
        }
        .map { bci_rds, bci_gene, ge_rds, _ge_gene -> 
            tuple(bci_rds, ge_rds, bci_gene) 
        }

    // RUN comparison BCI vs GE
    compareBCIwithGEvar(comparison_input_ch)

    // RUN metadata comparison plots (always runs, uses stub in personal_pc)
    // Prepare input for metadata plotting: combine comparison data, BCI clean data, and participant metadata
    metadata_plot_input_ch = compareBCIwithGEvar.out.BciVsGe_rds
        .join(cleanFormatBciPatientsVar.out.clean_rds, by: 0)  // Join by gene_name
        .combine(part_metadata_ch)
        .filter { gene_name, _comparison_rds, _bci_rds, part_metadata -> 
            part_metadata.name.contains(gene_name) 
        }
        .map { gene_name, comparison_rds, bci_rds, part_metadata -> 
            tuple(gene_name, comparison_rds, bci_rds, part_metadata) 
        }

    plotMetadataBCIvsGE(metadata_plot_input_ch)

    // RUN plotting QC results
    plotQCBciPatientsVar(cleanFormatBciPatientsVar.out.clean_rds, prot_files_ch, exon_files_ch)
            // out.plots: path(plot_file)

    // Collect versions
    ch_versions = ch_versions.mix(cleanFormatBciPatientsVar.out.versions)
    ch_versions = ch_versions.mix(compareBCIwithGEvar.out.versions)
    ch_versions = ch_versions.mix(plotMetadataBCIvsGE.out.versions) 
    ch_versions = ch_versions.mix(plotQCBciPatientsVar.out.versions)

    emit:
    clean_tsv           = cleanFormatBciPatientsVar.out.clean_tsv     // File TSV
    clean_rds           = cleanFormatBciPatientsVar.out.clean_rds     // path(rds_file)
    plots               = plotQCBciPatientsVar.out.plots              // PDFs
    comparison_tsv      = compareBCIwithGEvar.out.BciVsGe_tsv         // File TSV
    comparison_rds      = compareBCIwithGEvar.out.BciVsGe_rds         // path(rds_file)
    metadata_plots      = plotMetadataBCIvsGE.out.plots               // Metadata comparison PDFs
    versions            = ch_versions
}
