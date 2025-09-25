// Nextflow process to call the R script for filtering structural variants

process filterStructVar {
    tag { gene_name }
    label 'r_process'
    
    publishDir "${params.results_dir}/${gene_name}", mode: 'copy'

    input:
        tuple val(gene_name), path(clean_structvar_file), path(part_metadata_file)

    output:
        path("${gene_name}_struct_variants_filtered.tsv"), emit: filtered_clean_tsv
        tuple path("${gene_name}_struct_variants_filtered.rds"), val(gene_name), emit: filtered_clean_rds
        path("${gene_name}_struct_variants_filtered_stats.csv"), emit: stats_csv
        path "versions.yml", emit: versions

    script:
    """
    echo "Processing gene: ${gene_name} with annotation file: ${clean_structvar_file}"
    filterStructVar.R ${clean_structvar_file} ${gene_name} ${part_metadata_file}

    # Check output files created
    if [ ! -f "${gene_name}_struct_variants_filtered.tsv" ]; then
        echo "Error: TSV output file was not created"
        exit 1
    fi

    if [ ! -f "${gene_name}_struct_variants_filtered.rds" ]; then
        echo "Error: RDS output file was not created"
        exit 1
    fi

    if [ ! -f "${gene_name}_struct_variants_filtered_stats.csv" ]; then
        echo "Error: Stats output file was not created"
        exit 1
    fi
    
    echo "Both TSV and RDS files created successfully"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
    END_VERSIONS
    """
}