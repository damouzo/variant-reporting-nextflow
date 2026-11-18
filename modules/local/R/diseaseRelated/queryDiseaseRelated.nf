// Nextflow process for disease-related analysis (100kGP phenotype / ecds / hes_apc)

process queryDiseaseRelated {
    tag "${gene_name}_${filter_type}"
    label 'r_process_sql'

    publishDir "${params.results_dir}/${gene_name}/${filter_type}/disease_related", mode: 'copy'

    input:
        tuple val(gene_name), val(filter_type), path(partID_txt)
        val labkey_main
        val disease_terms
        val exclude_chapters

    output:
        path "${gene_name}_${filter_type}_disease_related_stats.csv", emit: stats_csv
        path "${gene_name}_${filter_type}_disease_related_participants.tsv", emit: participants_tsv
        path "${gene_name}_${filter_type}_disease_related_plots.pdf", emit: plots
        path "${gene_name}_${filter_type}_disease_related_comorbidity_*.pdf", emit: comorbidity_plots, optional: true
        path "versions.yml", emit: versions

    when:
        params.disease_related && params.disease_related.size() > 0 && params.enable_sql_queries

    script:
    def terms_arg = disease_terms.collect { it instanceof Map ? it.term : it }.join(',')
    def exclude_arg = exclude_chapters ? exclude_chapters.join(',') : 'Z,U,V,W,X,Y'
    """
    echo "Disease-related analysis: ${gene_name} / ${filter_type} -> ${terms_arg}"
    queryDiseaseRelated.R ${partID_txt} ${gene_name} ${filter_type} ${labkey_main} "${terms_arg}" "${exclude_arg}"

    if [ ! -f "${gene_name}_${filter_type}_disease_related_stats.csv" ]; then
        echo "Error: disease-related stats file was not created"
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    echo "STUB: disease-related analysis for ${gene_name} ${filter_type}"
    printf "term,ontology,n_filtered_total,n_filtered_with,pct_filtered,n_programme_total,n_programme_with,pct_programme,odds_ratio,ci_low,ci_high,p_value,fdr\\n" > ${gene_name}_${filter_type}_disease_related_stats.csv
    printf "participant_id\\tterm\\tin_filtered\\n" > ${gene_name}_${filter_type}_disease_related_participants.tsv
    echo "" > ${gene_name}_${filter_type}_disease_related_plots.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo "4.3.3")
    END_VERSIONS
    """
}
