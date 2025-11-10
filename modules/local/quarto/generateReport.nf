process generateReport {
    tag "$gene_name"
    label 'quarto_process'
    
    publishDir "${params.results_dir}/reports", mode: 'copy'

    input:
        tuple val(gene_name), val(report_type), path(template_file), path(input_files)

    output:
        path "${gene_name}_${report_type}.html", emit: html_report
        path "${gene_name}_${report_type}_files/", emit: report_files, optional: true
        path "versions.yml", emit: versions

    when:
        params.enable_quarto_report

    script:
    def input_files_list = input_files instanceof List ? input_files.join(' ') : input_files
    def report_name = "${gene_name}_${report_type}"
    """
    # Configure Quarto to use local writable directories
    export QUARTO_CACHE_DIR="\${PWD}/quarto_cache"
    export DENO_DIR="\${PWD}/deno_cache"
    export XDG_CACHE_HOME="\${PWD}/xdg_cache"

    # Create cache directories
    mkdir -p quarto_cache deno_cache xdg_cache data

    echo "Generating report: ${report_name}"
    echo "Using cache directories in: \${PWD}"

    # Copy template and replace placeholders with actual values
    cp ${template_file} ${report_name}.qmd

    # Replace placeholders - CUIDADO CON LAS COMILLAS
    sed -i 's/{{PROJECT_CODE}}/${params.project_code}/g' ${report_name}.qmd
    sed -i 's/{{PROJECT_NAME}}/${params.project_name.replaceAll("'", "\\\\'")}/g' ${report_name}.qmd
    sed -i 's/{{PROJECT_AUTHOR}}/${params.project_author}/g' ${report_name}.qmd
    sed -i 's/{{GENE_NAME}}/${gene_name}/g' ${report_name}.qmd

    # Debug: Show what was replaced
    echo "Replaced placeholders:"
    echo "PROJECT_CODE: ${params.project_code}"
    echo "PROJECT_NAME: ${params.project_name}"
    echo "PROJECT_AUTHOR: ${params.project_author}"
    echo "GENE_NAME: ${gene_name}"

    # Copy input files
    for file in ${input_files_list}; do
        if [ -f "\$file" ]; then
            cp "\$file" data/
            echo "Copied \$file to data/"
        fi
    done

    # List available files
    echo "Available data files:"
    ls -la data/

    # Render HTML
    echo "Rendering HTML report (PDF generation disabled for container optimization)..."
    quarto render ${report_name}.qmd --to html

    # Verify HTML output
    if [ -f "${report_name}.html" ]; then
        echo "HTML report generated successfully"
        echo "Note: PDF generation disabled to reduce container size from ~8GB to ~5GB"
    else
        echo "Error: HTML report was not created"
        exit 1
    fi

        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
    END_VERSIONS
    """
}