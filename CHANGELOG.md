# Changelog

## [v1.4.0] - 2026-11-09
- Per-filter participant outputs: `{gene}_small_variants_{filter}_participantID.txt` and `{gene}_small_variants_{filter}_participants_variants.tsv` (participants x variants with metadata)
- Disease-related analysis process (`params.disease_related`): queries 100kGP HPO phenotype, ecds (SNOMED) and hes_apc (ICD-10), computes prevalence and enrichment vs the programme, writes stats/plots under `{gene}/{filter}/disease_related/`
- Integrated both into the existing Quarto report

## [v1.3.0] - 2025-11-03
- Automated Quarto reporting system: Per-gene HTML reports with integrated variant statistics, QC plots, and filtering metrics
- Pipeline execution reports: Timeline visualization, resource usage tracking, DAG generation, and detailed trace files
- Advanced filtering framework: Gene-specific custom filtering strategies with configurable parameters

## [v1.2.0] - 2025-09-05
- Added new modules and workflows for BCI patients variants analysis, including data cleaning, QC plotting, and participant ID extraction.

## [v1.1.0] - 2025-08-29
- Added support for structural variant analysis, including new modules and workflows.
- Integrated HPC compatibility and Singularity container support for streamlined execution in diverse environments.
- Changed configuration system from YAML to custom config for improved flexibility.

## [v1.0.0] - 2025-08-25
- Finalized first version of all modules and subworkflows for small variant QC.
- Comprehensive documentation added, including pipeline execution reports and DAG visualization.
- Optimized R scripts for data cleaning and QC plotting.

## [v0.1.0] - 2025-08-20
- Enhanced `cleanFormatSmallVar.R` with additional data cleaning and output validation.
- Improved `plotQCsmallVar.R` to include advanced QC plots and error handling.
- Updated `main.nf` to include detailed logging and parameterized configurations.
- Refactored `nextflow.config` to dynamically load YAML configurations.
- Added `config_local.yaml` and `config_ge_local.yaml` for environment-specific settings.
- Refactored `cleanFormatSmallVar.nf` and `qc_small_variants.nf` to improve modularity and output handling.

## [Previous Changes]
- Initial project structure and documentation.
- Added Nextflow pipeline and R scripts for variant analysis.
- Included example data and results folders.
- Added README and .gitignore files.
