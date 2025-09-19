// Nextflow process to call the R script for plotting metadata comparison between BCI and GE patients

process plotMetadataBCIvsGE {
    tag { gene_name }
    label 'r_process'

    publishDir "${params.results_dir}/${gene_name}/plots/bci_patients_variants/comparison_BCI_GE", mode: 'copy'

    input:
    tuple val(gene_name), path(comparison_rds), path(bci_clean_rds), path(part_metadata_file)

    output:
    path "${gene_name}_*.pdf", emit: plots
    path "${gene_name}_*.csv", emit: csv    
    path "versions.yml", emit: versions

    script:
    """
    echo "Plotting metadata comparison for gene: ${gene_name}"
    echo "Input files:"
    echo "  Comparison data: ${comparison_rds}"
    echo "  BCI clean data: ${bci_clean_rds}"
    echo "  Participant metadata: ${part_metadata_file}"

    plotMetadataBCIvsGE.R ${gene_name} ${comparison_rds} ${bci_clean_rds} ${part_metadata_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    echo "STUB MODE: Creating mock files for ${gene_name}"

    # Create a mock CSV file
    echo "mock example when stub version" > ${gene_name}_mock_csv.csv

    # Create a simple mock PDF
    touch ${gene_name}_metadata_mockplot.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: "4.3.3"
    END_VERSIONS
    """
}