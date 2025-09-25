// Nextflow process to call the R script for extracting participant IDs from structural variants

process extractStructVarPartID {
    tag { gene_name }
    label 'r_process'

    publishDir "${params.results_dir}/${gene_name}", mode: 'copy'

    input:
        tuple val(gene_name), path(structvar_annot_file)
        val labkey_main

    output:
        tuple val(gene_name), path("${gene_name}_structural_variants_participantMetadata.rds"), emit: partMet_rds
        path("${gene_name}_structural_variants_participantID.txt"), emit: partID_txt
        path("${gene_name}_structural_variants_participantMetadata.tsv"), emit: partMet_tsv
        path "versions.yml", emit: versions

    script:
    """
    echo "Processing: ${gene_name} with annotation file: ${structvar_annot_file} and labkey: ${labkey_main}"
    extractStructVarPartID.R ${structvar_annot_file} ${gene_name} ${labkey_main}

    # Check output files created
    if [ ! -f "${gene_name}_structural_variants_participantID.txt" ]; then
        echo "Error: Participants ID output file was not created"
        exit 1
    fi

    if [ ! -f "${gene_name}_structural_variants_participantMetadata.tsv" ]; then
        echo "Error: Participants metadata output file was not created"
        exit 1
    fi

    echo "Participants SQL files created successfully"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    echo "STUB MODE: Using mock data for ${gene_name}"
    cp "${moduleDir}/mock_data/${gene_name}_structural_variants_participantID.txt" "${gene_name}_structural_variants_participantID.txt"     
    cp "${moduleDir}/mock_data/${gene_name}_structural_variants_participantMetadata.tsv" "${gene_name}_structural_variants_participantMetadata.tsv"
    echo "Mock files ready for ${gene_name}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo "4.3.3")
    END_VERSIONS
    """
}