// Nextflow script to call the R script for plotting filtered structural variants

process plotFilteredStructVar {
    tag { gene_name }
    label 'r_process'

    // Publicar todos los PDFs generados en subcarpeta organizada por gen
    publishDir "${params.results_dir}/${gene_name}/plots/structural_variants", mode: 'copy'

    input:
    tuple val(gene_name), path(filtered_table), path(prot_file), path(exon_file), path(part_metadata_file)

    // Save all PDFs generated
    output:
    path "*.pdf", emit: filtered_plots
    path "versions.yml", emit: versions

    script:
    """
    echo "Plotting filtered variants for gene: ${gene_name}"
    plotFilteredStructVar.R ${filtered_table} ${gene_name} ${prot_file} ${exon_file} ${part_metadata_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
    END_VERSIONS
    """
}