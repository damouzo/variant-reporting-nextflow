// Nextflow script to call the R script for plotting QC of small variants

process plotQCsmallVar {
    tag { gene_name }
    label 'r_process'

    // Publicar todos los PDFs generados en subcarpeta organizada por gen
    publishDir "${params.results_dir}/${gene_name}/plots/small_variants", mode: 'copy'

    input:
    tuple path(clean_table), val(gene_name), path(prot_file), path(exon_file)

    // Save all PDFs generated
    output:
    path "*.pdf"
    path "versions.yml", emit: versions

    script:
    """
    echo "Plotting QC for gene: ${gene_name}"
    plotQCsmallVar.R ${clean_table} ${gene_name} ${prot_file} ${exon_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
    END_VERSIONS
    """
}