// Nextflow process to call the R script for filtering small variants

process filterSmallVar {
    tag { gene_name }
    label 'r_process'
    
    publishDir "${params.results_dir}/${gene_name}", mode: 'copy'

    input:
        tuple path(clean_smallvar_file), val(gene_name)

    output:
        path("${gene_name}small_variants_filtered.tsv"), emit: filtered_clean_tsv
        tuple path("${gene_name}small_variants_filtered.rds"), val(gene_name), emit: filtered_clean_rds
        path "versions.yml", emit: versions

    script:
    """
    echo "Processing gene: ${gene_name} with annotation file: ${clean_smallvar_file}"
    filterSmallVar.R ${clean_smallvar_file} ${gene_name}

    # Check output files created
    if [ ! -f "${gene_name}small_variants_filtered.tsv" ]; then
        echo "Error: TSV output file was not created"
        exit 1
    fi

    if [ ! -f "${gene_name}small_variants_filtered.rds" ]; then
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