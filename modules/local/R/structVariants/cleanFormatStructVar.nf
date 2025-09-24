// Nextflow proces to call the R script for cleaning and formatting structural variants

process cleanFormatStructVar {
    tag { gene_name }
    label 'r_process'
    
    publishDir "${params.results_dir}/${gene_name}", mode: 'copy'

    input:
        tuple path(struct_cnv_variants), path(struct_sv_variants), val(gene_name)

    output:
        tuple path("${gene_name}_structural_variants.rds"), val(gene_name), emit: clean_rds
        path("${gene_name}_structural_variants.tsv")
        path "versions.yml", emit: versions

    script:
    """
    echo "Processing gene: ${gene_name} with annotation files:"
    echo "  CNV: ${struct_cnv_variants}"
    echo "  SV:  ${struct_sv_variants}"
    cleanFormatStructVar.R ${struct_cnv_variants} ${struct_sv_variants} ${gene_name}
    
    # Check output files created
    if [ ! -f "${gene_name}_structural_variants.rds" ]; then
        echo "Error: Structural variants RDS output file was not created"
        exit 1
    fi

    if [ ! -f "${gene_name}_structural_variants.tsv" ]; then
        echo "Error: Structural variants TSV output file was not created"
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
    END_VERSIONS
    """
}