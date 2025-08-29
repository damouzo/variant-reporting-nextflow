# Changelog

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
