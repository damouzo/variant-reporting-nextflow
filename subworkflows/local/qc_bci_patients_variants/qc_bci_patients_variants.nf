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
    clean_variants_ch       // channel: GE small variants from SmallVar subworkflow
    part_metadata_ch        // channel: participant metadata from SmallVar subworkflow

    main:
    // Create versions channel
    ch_versions = Channel.empty()

    // Load reference annotation files
    prot_files_ch = Channel.fromPath("${params.prot_dir}/*.gff")
    exon_files_ch = Channel.fromPath("${params.exon_dir}/*.tsv")

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

        // Clean and format BCI patients variants (same process, different input)
        cleanFormatBciPatientsVar(pre_annotated_ch)
    }
    
    // Combine BCI clean data with clean variants from the first subworkflow
    comparison_input_ch = cleanFormatBciPatientsVar.out.clean_rds
        .join(clean_variants_ch)
        .map { gene_name, bci_rds, ge_rds -> 
            tuple(bci_rds, ge_rds, gene_name) 
        }

    // RUN comparison BCI vs GE
    compareBCIwithGEvar(comparison_input_ch)

    // Prepare input for metadata plotting
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

    metadata_plot_input_ch = compareBCIwithGEvar.out.BciVsGe_rds
    .join(cleanFormatBciPatientsVar.out.clean_rds)  // [gene_name, comparison_rds, bci_rds]
    .map { gene_name, comparison_rds, bci_rds -> 
        tuple(gene_name, comparison_rds, bci_rds)
    }
    .combine(part_metadata_ch)  // [gene_name, comparison_rds, bci_rds, metadata_file]
    .filter { gene_name, _comparison_rds, _bci_rds, metadata_file -> 
        metadata_file.name.contains(gene_name)
    }
    .map { gene_name, comparison_rds, bci_rds, metadata_file -> 
        tuple(gene_name, comparison_rds, bci_rds, metadata_file)
    }

    // RUN metadata comparison plots
    plotMetadataBCIvsGE(metadata_plot_input_ch)
        // out.plots: path(plot_file)
        // out.csv: path(csv_file)

    plot_qc_input_ch = cleanFormatBciPatientsVar.out.clean_rds
        .join(prot_files_with_gene_ch)
        .join(exon_files_with_gene_ch)
        .map { gene_name, clean_rds, prot_file, exon_file -> 
            tuple(gene_name, clean_rds, prot_file, exon_file) 
        }

    // RUN plotting QC results
    plotQCBciPatientsVar(plot_qc_input_ch)
            // out.plots: path(plot_file)

    // Collect versions
    ch_versions = ch_versions.mix(cleanFormatBciPatientsVar.out.versions)
    ch_versions = ch_versions.mix(compareBCIwithGEvar.out.versions)
    ch_versions = ch_versions.mix(plotMetadataBCIvsGE.out.versions) 
    ch_versions = ch_versions.mix(plotQCBciPatientsVar.out.versions)

    emit:
    clean_tsv           = cleanFormatBciPatientsVar.out.clean_tsv     // File TSV
    clean_rds           = cleanFormatBciPatientsVar.out.clean_rds     // tuple(val(gene_name), path(rds_file))
    plots               = plotQCBciPatientsVar.out.plots              // PDFs
    comparison_tsv      = compareBCIwithGEvar.out.BciVsGe_tsv         // File TSV
    comparison_rds      = compareBCIwithGEvar.out.BciVsGe_rds         // tuple(val(gene_name), path(rds_file))
    metadata_plots      = plotMetadataBCIvsGE.out.plots               // Metadata comparison PDFs
    versions            = ch_versions
}
