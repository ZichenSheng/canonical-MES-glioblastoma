# Glioblastoma manuscript reproducibility

This repository provides the minimal analytical code and aggregate validation supporting the principal quantitative results of **Transcriptomic signatures preserve biological direction without uniquely specifying state in glioblastoma**. Raw study datasets are not redistributed and must be obtained from their original sources. This is a bounded manuscript reproducibility release, not a one-command reconstruction from raw data.

The current Results authority is the September 19, 2026 manuscript, with 78 claim-level entries and six main figures. The historical v1.0.0 is preserved; v2 follows the current six-section manuscript. See [authority](docs/CURRENT_RESULTS_AUTHORITY.md), [claim crosswalk](docs/RESULT_CODE_CROSSWALK.tsv) and [figure map](docs/RESULT_MAP.md).

| Module | Actual certified capability |
|---|---|
| 01_inspect_signatures.R | Inspect MES95 and Primary201 membership definitions |
| 02_patient_external_transfer.py | Recompute patient-exclusive ecology transfer and patient-bootstrap CI from external prepared arrays |
| 03_anatomical_resolution.py | Recompute tumour-blocked association and leave-one-tumour-out resolution from external prepared scores |
| 04_longitudinal_deformation.py | Recompute patient-exclusive all-k deformation from external paired scores |
| 05_design_capability.R | Validate exact-support qualification and recovery from frozen aggregate counts |
| 06_crc04_exact_inference.R | Recompute frozen measurement and exact inference from downloaded public processed WTA RDS; includes synthetic fixture |
| 07_sampling_composition.py | Recompute CARE sampling-null comparisons and composition decomposition from external frozen prepared inputs |
| 08_aggregate_evidence.py | Recalculate overlap, localization, state-resolution, stability and recovery summaries from safe aggregate inputs |

Eight principal scripts have PASS receipts for their stated scope. The crosswalk distinguishes aggregate validation, public prepared-input reproduction and source-dependent execution; a numerical PASS does not certify omitted preprocessing. No claim is labelled FULL_FROM_SOURCE_DATA. Graphical assembly and real-data plot reconstruction are outside this release.

Quick start (R 4.4.1; Python 3.11):

```sh
Rscript -e 'renv::restore(prompt=FALSE)'
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-lock.txt
MANUSCRIPT_PYTHON=.venv/bin/python Rscript analysis/99_validate_release.R
```

Install renv separately if needed. The validator runs the immediately runnable modules and checks frozen table hashes, schemas, printed values, numerical receipts, denominator and estimand locks. It exits nonzero on mismatch. [Data access](docs/DATA_ACCESS.md) documents Tier 2 CRC04 and Tier 3 source-dependent modes; [datasets.tsv](config/datasets.tsv) registers 18 current-paper resource/subset entries.

R1 uses 57 ALL_AUTHOR_ANNOTATED regional measurements from 12 patients. R4 uses 57 MALIGNANT_ONLY regional pseudobulks from 12 donors. The approved [Methods wording correction](docs/MANUSCRIPT_DISCREPANCY_LOG.tsv) changes no result or executed analysis. CRC04's historical four intermediate-vector exact-equality discrepancies remain disclosed; its exact inferential chain reproduces.

The independent [evaluably 0.9.0 release](https://github.com/ZichenSheng/evaluably/releases/tag/v0.9.0), DOI [10.5281/zenodo.22846293](https://doi.org/10.5281/zenodo.22846293), is frozen and is not vendored or an execution dependency here. That DOI is not the manuscript-code DOI. The v2 version DOI will be added to main after archival without moving the release tag.
