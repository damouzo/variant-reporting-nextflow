// Nextflow process to call the R script for cleaning and formatting BCI patients variants

process compareBCIwithGEvar {
    tag { gene_name }
    label 'r_process'
    
    publishDir "${params.results_dir}/${gene_name}", mode: 'copy'

    input:
        tuple path(bci_SmallVar_rds), path(ge_SmallVar_rds), val(gene_name)

    output:
        tuple val(gene_name), path("${gene_name}_compared_bci_ge_small_variants.tsv"), emit: BciVsGe_tsv
        tuple val(gene_name), path("${gene_name}_compared_bci_ge_small_variants.rds"), emit: BciVsGe_rds
        path "versions.yml", emit: versions

    script:
    """
    echo "Processing: ${gene_name} with BCI file: ${bci_SmallVar_rds}, and GE file: ${ge_SmallVar_rds}"
    compareBCIwithGEvar.R ${bci_SmallVar_rds} ${ge_SmallVar_rds} ${gene_name}

    # Check output files created
    if [ ! -f "${gene_name}_compared_bci_ge_small_variants.tsv" ]; then
        echo "Error: TSV output file was not created"
        exit 1
    fi

    if [ ! -f "${gene_name}_compared_bci_ge_small_variants.rds" ]; then
        echo "Error: RDS output file was not created"
        exit 1
    fi
    
    echo "Both TSV and RDS files created successfully"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
    END_VERSIONS
    """
}