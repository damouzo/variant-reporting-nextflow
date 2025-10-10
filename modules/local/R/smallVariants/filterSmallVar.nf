// Nextflow process to call the R script for filtering small variants

process filterSmallVar {
    tag { "${gene_name}_${filter_type}" }
    label 'r_process'
    
    publishDir "${params.results_dir}/${gene_name}/${filter_type}", mode: 'copy'

    input:
        tuple val(gene_name), path(clean_smallvar_file), path(part_metadata_file), val(filter_type), val(filter_config)

    output:
        path "${gene_name}_small_variants_filtered_stats.csv", emit: stats_csv
        path "${gene_name}_small_variants_filtered.tsv", emit: filtered_clean_tsv
        tuple val(gene_name), val(filter_type), path("${gene_name}_small_variants_filtered.rds"), emit: filtered_clean_rds
        path "versions.yml", emit: versions

    script:
    def filter_args = filter_config ? filter_config.collect { key, value -> "--${key} ${value}" }.join(' ') : ''
    """
    echo "Processing gene: ${gene_name} with filter type: ${filter_type}"
    echo "Filter configuration: ${filter_config}"
    
    filterSmallVar.R ${clean_smallvar_file} ${gene_name} ${part_metadata_file} ${filter_type} ${filter_args}

    # Check output files created
    if [ ! -f "${gene_name}_small_variants_filtered.tsv" ]; then
        echo "Error: TSV output file was not created"
        exit 1
    fi

    if [ ! -f "${gene_name}_small_variants_filtered.rds" ]; then
        echo "Error: RDS output file was not created"
        exit 1
    fi

    if [ ! -f "${gene_name}_small_variants_filtered_stats.csv" ]; then
        echo "Error: Stats output file was not created"
        exit 1
    fi
    
    echo "Both TSV and RDS files created successfully"
    echo "Filtered files created successfully for ${filter_type}"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
    END_VERSIONS
    """
}