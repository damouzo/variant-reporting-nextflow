// Nextflow script to call the R script for plotting QC of BCI patients variants

process plotQCBciPatientsVar {
    tag { gene_name }
    label 'r_process'

    // Publicar todos los PDFs generados en subcarpeta organizada por gen
    publishDir "${params.results_dir}/${gene_name}/plots/bci_patients_variants", mode: 'copy'

    input:
    tuple val(gene_name), path(bci_rds), path(prot_file), path(exon_file)

    // Save all PDFs generated
    output:
    path "*.pdf", emit: plots
    path "versions.yml", emit: versions

    script:
    """
    echo "Plotting BCI patients variants QC for gene: ${gene_name}"
    plotQCBciPatientsVar.R ${bci_rds} ${gene_name} ${prot_file} ${exon_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
    END_VERSIONS
    """
}