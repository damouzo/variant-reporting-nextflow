// normalizeCombineTSV.nf

process normalizeCombineTSV {
    tag "${gene_name}"
    label 'r_process'

    publishDir "${params.bci_annotated_dir}", mode: 'copy', pattern: "*_normalized_combined_data.tsv"

    input:
    tuple val(gene_name), path(bci_patients_var_files)

    output:
    tuple val(gene_name), path("${gene_name}_normalized_combined_data.tsv"), emit: normalized_tsv
    tuple val(gene_name), path("${gene_name}_vep_input.txt"), emit: vep_input
    path "versions.yml", emit: versions

    script:
    """
    echo "Processing gene: ${gene_name} with files: ${bci_patients_var_files.join(', ')}"
    normalizeCombineTSV.R ${gene_name} ${bci_patients_var_files.join(' ')}

    # Check output file created
    if [ ! -f "${gene_name}_normalized_combined_data.tsv" ]; then
        echo "Error: Normalized combined TSV file was not created"
        exit 1
    fi

    echo "Normalized and combined TSV file created successfully"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
    END_VERSIONS
    """
}