// Nextflow proces to call the R script for cleaning and formatting structural variants

process cleanFormatStructVar {
    tag { gene_name }
    //conda "r-base=4.3.0"
    publishDir "${params.results_dir}/${gene_name}", mode: 'copy'

    input:
        tuple path(struct_cnv_variants), path(struct_sv_variants), val(gene_name)

    output:
        tuple path("${gene_name}_structural_variants.rds"), val(gene_name), emit: struct_rds
        val("${gene_name}_structural_variants.tsv")
        path "versions.yml", emit: versions

    script:
    """
    echo "Processing gene: ${gene_name} with annotation files:"
    echo "  CNV: ${struct_cnv_variants}"
    echo "  SV:  ${struct_sv_variants}"
    cleanFormatStructVar.R ${struct_cnv_variants} ${struct_sv_variants} ${gene_name}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch ${gene_name}_structural_variants.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.0
    END_VERSIONS
    """
}