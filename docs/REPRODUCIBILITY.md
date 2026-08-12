# Reproducibility guide

## Scope

This repository preserves the final analytical code supporting the manuscript's principal conclusions. It is a minimal code release, not a complete project-history archive or a bundled data distribution. The statistical calculations are retained; local I/O paths and internal run-management language were replaced for public release.

## Configuration

Set three environment variables from the repository root:

```bash
export GBM_MES_DATA_ROOT=/path/to/local/data
export GBM_MES_OUTPUT_ROOT=results
export GBM_MES_REPO_ROOT="$PWD"
```

`config/paths.example.yaml` documents the same roots for users who manage paths with YAML. Scripts use the environment variables directly.

## Prepared-input layout

The code refers to public datasets under `GBM_MES_DATA_ROOT`, grouped into `bulk/`, `single_cell/`, `spatial/`, `clinical/`, `depmap/`, and `prepared/`. The `prepared/` directory contains small cohort mappings, verified sample annotations, and interoperable tables derived from the public downloads. These patient-level inputs are deliberately not included in the repository.

Expected filenames are explicit near the top of each script. This makes the I/O contract inspectable without embedding a machine-specific path. Users adapting a public download may either use those filenames or update only the corresponding I/O declarations.

## Analysis order

1. `scripts/01_scoring`: canonical and comparator scoring in pseudo-bulk and bulk cohorts.
2. `scripts/02_single_cell`: patient-level malignant and compartment analyses.
3. `scripts/03_decomposition`: BayesPrism reuse, state-specific controls, and bulk–malignant discordance.
4. `scripts/04_patient_realization`: patient fingerprints, residual structure, and cross-cohort rank-2 subspace validation.
5. `scripts/05_provenance`: gene-origin contributions, ablations, stability, ICB transport, and DepMap controls.
6. `scripts/06_spatial`: Visium and Ivy GAP context analyses plus the formal spatial meta-test.
7. `scripts/07_clinical`: ICB, longitudinal, GLASS, non-proportional-hazards, and sensitivity analyses.

Modules reuse outputs from earlier modules where their paths point into `GBM_MES_OUTPUT_ROOT`. Full-data execution can be computationally intensive and was not repeated during public-release preparation. Python compilation and R parse checks are the release-level static verification.

## Figures

The retained scripts generate analytical tables and panel-level quantities. Final figure assembly and graphical annotation were performed separately. The repository therefore does not claim pixel-identical reconstruction of the submitted figures.
