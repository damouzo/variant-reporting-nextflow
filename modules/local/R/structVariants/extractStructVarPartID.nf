// Nextflow process to call the R script for extracting participant IDs from structural variants

process extractStructVarPartID {
    tag { gene_name }
    label 'r_process'
    label 'sql_access'

    publishDir "${params.results_dir}/${gene_name}", mode: 'copy'

    when: 
    params.enable_sql_queries

    input:
        tuple path(strucvar_annot_file), val(gene_name)
        val labkey_main

    output:
        path("${gene_name}_structural_variants_participantID.txt"), emit: partID
        path("${gene_name}_structural_variants_participantMetadata.tsv"), emit: partMet
        path "versions.yml", emit: versions

    script:
    """
    echo "Processing: ${gene_name} with annotation file: ${strucvar_annot_file} and labkey: ${labkey_main}"
    extractStructVarPartID.R ${strucvar_annot_file} ${gene_name} ${labkey_main}

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
}