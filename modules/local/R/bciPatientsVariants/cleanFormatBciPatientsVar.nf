// Nextflow process to call the R script for cleaning and formatting BCI patients variants

process cleanFormatBciPatientsVar {
    tag { gene_name }
    label 'r_process'
    
    publishDir "${params.results_dir}/${gene_name}", mode: 'copy'

    input:
        tuple val(gene_name), path(vep_output), path(original_data)

    output:
        path("${gene_name}_bci_patients_variants.tsv"), emit: clean_tsv
        tuple path("${gene_name}_bci_patients_variants.rds"), val(gene_name), emit: clean_rds
        path "versions.yml", emit: versions

    script:
    """
    echo "Processing: ${gene_name} with BCI patients variants file: ${vep_output}"
    cleanFormatBciPatientsVar.R ${gene_name} ${vep_output} ${original_data}

    # Check output files created
    if [ ! -f "${gene_name}_bci_patients_variants.tsv" ]; then
        echo "Error: TSV output file was not created"
        exit 1
    fi

    if [ ! -f "${gene_name}_bci_patients_variants.rds" ]; then
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