// Nextflow process to call the R script for cleaning and formatting small variants

process cleanFormatSmallVar {
    tag { gene_name }
    label 'r_process'
    
    publishDir "${params.results_dir}/${gene_name}", mode: 'copy'

    input:
        tuple path(smallvar_annot_file), val(gene_name)

    output:
        path("${gene_name}_small_variants.tsv"), emit: clean_tsv
        tuple path("${gene_name}_small_variants.rds"), val(gene_name), emit: clean_rds
        path "versions.yml", emit: versions

    script:
    """
    echo "Processing gene: ${gene_name} with annotation file: ${smallvar_annot_file}"
    cleanFormatSmallVar.R ${smallvar_annot_file} ${gene_name}

    # Check output files created
    if [ ! -f "${gene_name}_small_variants.tsv" ]; then
        echo "Error: TSV output file was not created"
        exit 1
    fi
    
    if [ ! -f "${gene_name}_small_variants.rds" ]; then
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