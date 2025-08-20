// Nextflow script to call the R script for plotting QC of structural variants

process plotQCstructVar {
    tag { gene_name }

    // Publicar todos los PDFs generados en subcarpeta organizada por gen
    publishDir "${params.results_dir}/${gene_name}/plots/struct_variants", mode: 'copy'

    input:
    tuple path(clean_table), val(gene_name)

    // Save all PDFs generated
    output:
    path "plots/struct_variants/*.pdf"

    script:
    """
    echo "Plotting QC for gene: ${gene_name}"
    plotQCstructVar.R ${clean_table} ${gene_name}
    """
}