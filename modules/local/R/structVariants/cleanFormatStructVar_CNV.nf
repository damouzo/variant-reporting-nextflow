// Nextflow process to call the R script for cleaning and formatting structural CNV variants

process cleanFormatStructVar_CNV {
    tag { gene_name }
    label 'r_process'

    input:
        tuple path(cnv_germline_annot_file), path(cnv_somatic_annot_file), val(gene_name)

    output:
        tuple path("${gene_name}_struct_cnv_variants.rds"), val(gene_name), emit: clean_cnv
        path "versions.yml", emit: versions

    script:
    """
    echo "Processing gene: ${gene_name}"
    echo "Input files:"
    echo "  Germline: ${cnv_germline_annot_file}"
    echo "  Somatic:  ${cnv_somatic_annot_file}"
    
    # Run the R script
    cleanFormatStructVar_CNV.R ${cnv_germline_annot_file} ${cnv_somatic_annot_file} ${gene_name}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
    END_VERSIONS
    """
}