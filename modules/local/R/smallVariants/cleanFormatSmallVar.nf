// Nextflow proces to call the R script for cleaning and formatting small variants

process cleanFormatSmallVar {
    tag { gene_name }
    publishDir "${params.results_dir}/${gene_name}", mode: 'copy'

    input:
        tuple path(smallvar_annot_file), val(gene_name)

    output:
        tuple path("${gene_name}_small_variants.tsv"), val(gene_name)

    script:
    """
    echo "Processing gene: ${gene_name} with annotation file: ${smallvar_annot_file}"
    CleanFormatSmallVar.R ${smallvar_annot_file} ${gene_name}
    """
}