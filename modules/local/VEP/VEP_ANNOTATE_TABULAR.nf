process VEP_ANNOTATE_TABULAR {
    tag "${gene_name}"
    label 'vep_process'

    input:
    tuple val(gene_name), path(vep_input)
    path vep_cache

    output:
    tuple val(gene_name), path("${gene_name}_vep_output.txt"), emit: vep_output
    path "${gene_name}_vep_stats.html", emit: report
    path "versions.yml", emit: versions

    script:
    """
    vep \
        --input_file ${vep_input} \
        --output_file ${gene_name}_vep_output.txt \
        --species homo_sapiens \
        --assembly GRCh38 \
        --cache_dir ${vep_cache} \
        --offline \
        --format hgvs \
        --tab \
        --everything \
        --shift_hgvs 1 \
        --check_existing \
        --stats_file ${gene_name}_vep_stats.html \
        --fork ${task.cpus} \
        --force_overwrite

    # Verificar que se generó el output
    if [ ! -f "${gene_name}_vep_output.txt" ]; then
        echo "Error: VEP output file was not created"
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(\$VEP_CMD --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}