// Nextflow process to call the R script for plotting metadata comparison between BCI and GE patients

process plotMetadataBCIvsGE {
    tag { gene_name }
    label 'r_process'

    publishDir "${params.results_dir}/${gene_name}/plots/bci_patients_variants/comparison_BCI_GE", mode: 'copy'

    input:
    tuple val(gene_name), path(comparison_rds), path(bci_clean_rds), path(part_metadata_file)

    output:
    path "*.pdf", emit: plots
    path "versions.yml", emit: versions

    script:
    """
    echo "Plotting metadata comparison for gene: ${gene_name}"
    echo "Input files:"
    echo "  Comparison data: ${comparison_rds}"
    echo "  BCI clean data: ${bci_clean_rds}"
    echo "  Participant metadata: ${part_metadata_file}"
    
    plotMetadataBCIvsGE.R ${comparison_rds} ${bci_clean_rds} ${part_metadata_file} ${gene_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    echo "STUB MODE: Creating mock plot for ${gene_name}"
    
    # Create a simple mock PDF
    echo "%PDF-1.4" > ${gene_name}_metadata_mockplot.pdf
    echo "1 0 obj" >> ${gene_name}_metadata_mockplot.pdf
    echo "<<" >> ${gene_name}_metadata_mockplot.pdf
    echo "/Type /Catalog" >> ${gene_name}_metadata_mockplot.pdf
    echo "/Pages 2 0 R" >> ${gene_name}_metadata_mockplot.pdf
    echo ">>" >> ${gene_name}_metadata_mockplot.pdf
    echo "endobj" >> ${gene_name}_metadata_mockplot.pdf
    echo "Mock metadata plot created for ${gene_name}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo "4.3.3")
    END_VERSIONS
    """
}