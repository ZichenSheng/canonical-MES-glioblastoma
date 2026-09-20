# Canonical MES programs in glioblastoma

Code accompanying the manuscript on the biological interpretation, patient-level realization, and ecological context of canonical mesenchymal programs in glioblastoma.

## Overview

This repository contains the analysis code used for the principal results reported in the manuscript. The analyses evaluate canonical mesenchymal (MES) scoring across bulk, single-cell, spatial, decomposition, patient-realization, provenance, and clinical settings. The central interpretation is that canonical MES is informative but non-identifying and context-sensitive: the score captures reproducible biology without uniquely specifying one cellular source or one patient-level realization. Raw study datasets are not redistributed.

## Repository contents

- `scripts/`: 37 final analysis scripts organized by manuscript module.
- `resources/signatures/`: small, published or prespecified gene-set definitions used by the retained analyses.
- `config/`: example path configuration without machine-specific locations.
- `reproducibility/`: public-safe dataset registry and recorded software versions.
- `docs/`: data-access and reproducibility details.

The script-to-result map is summarized in [`scripts/README.md`](scripts/README.md).

| Manuscript figure | Main analytical scope |
|---|---|
| Figure 1 | Bulk MES ecology and benchmarking |
| Figure 2 | Cellular-source architecture and provenance |
| Figure 3 | MES-conditioned biological realization, CARE validation, and construct tests |
| Figure 4 | Gene provenance and DepMap model-context analysis |
| Figure 5 | Spatial, anatomical, and histopathological context |
| Figure 6 | Clinical context, specificity, and evidence boundaries |

## Data availability

Raw datasets are not redistributed in this repository. Public accessions and source pages are listed in [`reproducibility/dataset_registry.tsv`](reproducibility/dataset_registry.tsv) and described in [`docs/DATA_AVAILABILITY.md`](docs/DATA_AVAILABILITY.md). Each dataset remains subject to its original terms of use.

## Reproducibility

Set local roots before running a module:

```bash
export GBM_MES_DATA_ROOT=/path/to/local/data
export GBM_MES_OUTPUT_ROOT=results
export GBM_MES_REPO_ROOT="$PWD"
```

The repository records the final statistical code but does not claim a one-command reconstruction from raw data. Several modules require public datasets to be downloaded and organized into the prepared-input layout described in [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md). Final figure assembly and graphical annotation were performed separately from the analytical scripts.

## Software

The final analyses used R 4.4.1 and Python 3.11.8. Major package versions are recorded in [`reproducibility/software_versions.md`](reproducibility/software_versions.md); Python requirements are listed in [`requirements.txt`](requirements.txt).

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff).

## License

Code in this repository is released under the [MIT License](LICENSE). The license does not apply to third-party datasets, which remain subject to their original terms of use.

## Frozen manuscript-associated release

v1.0.0: [doi:10.5281/zenodo.21962295](https://doi.org/10.5281/zenodo.21962295)
