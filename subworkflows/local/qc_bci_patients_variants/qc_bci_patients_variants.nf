// qc_bci_patients_variants.nf

include { TSV_TO_VEP_INPUT } from '../../../modules/local/VEP/TSV_TO_VEP_INPUT.nf'
include { VEP_ANNOTATE_TABULAR } from '../../../modules/local/VEP/VEP_ANNOTATE_TABULAR.nf'
include { cleanFormatBciPatientsVar } from '../../../modules/local/R/bciPatientsVariants/cleanFormatBciPatientsVar.nf'
include { compareBCIwithGEvar } from '../../../modules/local/R/bciPatientsVariants/compareBCIwithGEvar.nf'
include { plotQCBciPatientsVar } from '../../../modules/local/R/bciPatientsVariants/plotQCBciPatientsVar.nf'


workflow QC_BCI_PATIENTS_VARIANTS {
    take:
    bci_patients_var_ch     // channel: BCI patients variant files
    gene_list_ch            // channel: list of genes
    prot_files_ch           // channel: protein annotation files
    exon_files_ch           // channel: exon annotation files
    clean_variants_ch       // channel: GE small variants from SmallVar subworkflow

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

        // Run VEP annotation
        VEP_ANNOTATE_TABULAR(TSV_TO_VEP_INPUT.out.vep_input)

        // Merge VEP annotations with original data
        merge_input_ch = VEP_ANNOTATE_TABULAR.out.vep_output
            .join(TSV_TO_VEP_INPUT.out.original_data)

        // Clean and format BCI patients variants
        cleanFormatBciPatientsVar(merge_input_ch)

        // Collect VEP versions
        ch_versions = ch_versions.mix(TSV_TO_VEP_INPUT.out.versions)
        ch_versions = ch_versions.mix(VEP_ANNOTATE_TABULAR.out.versions)

    } else {
        // Pre-annotated data path (for ge_pc)
        log.info "Using pre-annotated BCI patients variants from ${params.bci_annotated_dir}"
        
        // Use pre-annotated files
        pre_annotated_ch = Channel
            .fromPath("${params.bci_annotated_dir}/*.tsv")
            .combine(gene_list_ch)
            .filter { file, gene -> file.name.contains(gene) }
            .map { file, gene -> tuple(gene, file, file) // gene, vep_output, original_data
            }

        // Clean and format BCI patients variants (same process, different input)
        cleanFormatBciPatientsVar(pre_annotated_ch)
    }
    
    // Combine BCI clean data with clean variants from the first subworkflow
    comparison_input_ch = cleanFormatBciPatientsVar.out.clean_rds
        .combine(clean_variants_ch)
        .filter { bci_rds, bci_gene, ge_rds, ge_gene -> 
            bci_gene == ge_gene 
        }
        .map { bci_rds, bci_gene, ge_rds, ge_gene -> 
            tuple(bci_rds, ge_rds, bci_gene) 
        }

    // RUN comparison BCI vs GE
    compareBCIwithGEvar(comparison_input_ch)

    // RUN plotting QC results
    plotQCBciPatientsVar(cleanFormatBciPatientsVar.out.clean_rds, prot_files_ch, exon_files_ch)
            // out.plots: path(plot_file)

    // Collect versions
    ch_versions = ch_versions.mix(TSV_TO_VEP_INPUT.out.versions)
    ch_versions = ch_versions.mix(VEP_ANNOTATE_TABULAR.out.versions)
    ch_versions = ch_versions.mix(cleanFormatBciPatientsVar.out.versions)
    ch_versions = ch_versions.mix(compareBCIwithGEvar.out.versions)
    ch_versions = ch_versions.mix(plotQCBciPatientsVar.out.versions)

    emit:
    clean_tsv           = cleanFormatBciPatientsVar.out.clean_tsv     // File TSV
    clean_rds           = cleanFormatBciPatientsVar.out.clean_rds     // path(rds_file)
    plots               = plotQCBciPatientsVar.out.plots              // PDFs
    comparison_tsv      = compareBCIwithGEvar.out.BciVsGe_tsv         // File TSV
    comparison_rds      = compareBCIwithGEvar.out.BciVsGe_rds         // path(rds_file)
    versions            = ch_versions
}
