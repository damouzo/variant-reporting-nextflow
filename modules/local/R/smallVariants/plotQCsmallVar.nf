// Nextflow script to call the R script for plotting QC of small variants

process plotQCsmallVar {
    tag { gene_name }

    // Publicar todos los PDFs generados en subcarpeta organizada por gen
    publishDir "${params.results_dir}/${gene_name}/plots/small_variants", mode: 'copy'

    input:
    tuple path(clean_table), val(gene_name), path(prot_file), path(exon_file)

    // Save all PDFs generated
    output:
    path "plots/small_variants/*.pdf"

    script:
    """
    echo "Plotting QC for gene: ${gene_name}"
    plotQCsmallVar.R ${clean_table} ${gene_name} ${prot_file} ${exon_file}
    """
}