// Nextflow process to call the R script for cleaning and formatting structural SV variants

process cleanFormatStructVar_SV {
    tag { gene_name }
    //conda "r-base=4.3.0"

    input:
        tuple path(sv_germline_annot_file), path(sv_somatic_annot_file), val(gene_name)

    output:
        tuple path("${gene_name}_struct_sv_variants.rds"), val(gene_name)
        path "versions.yml", emit: versions

    script:
    """
    echo "Processing gene: ${gene_name}"
    echo "Input files:"
    echo "  Germline: ${sv_germline_annot_file}"
    echo "  Somatic:  ${sv_somatic_annot_file}"

    # Run the R script
    cleanFormatStructVar_SV.R ${sv_germline_annot_file} ${sv_somatic_annot_file} ${gene_name}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch ${gene_name}_struct_sv_variants.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.0
    END_VERSIONS
    """
}