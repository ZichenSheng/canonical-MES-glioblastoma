# Data access and prepared-input contracts

The current-paper registry is config/datasets.tsv. No study observation rows, expression matrices, donor capacities, simulated world-instance tables, images, or source archives are included. Aggregate tables are published reference outputs, not substitutes for raw-data regeneration.

Tier 1: scripts 01 and 05, the synthetic default of 06, and 99 --integrity-only run without source databases. Install the minimal R lock (digest 0.6.39; R 4.4.1).

Tier2 is certified for CRC04 using the source processed WTA RDS as described below. For the other analyses, obtaining the original source downloads does not by itself certify the exact prepared arrays; they remain Tier3.

Tier3: scripts02–04 and07 accept externally held prepared data. Numerical kernels were tested on the current frozen inputs locally, and only aggregate receipts are included.

* 02 input is an NPZ with scores (57 x 200), endpoints (57 x 4 ordered MES, HYPOXIA, MATRIX, MYELOID), groups (57 strings representing 12 patients), and model_ids (200). Original scoring uses the ALL_AUTHOR_ANNOTATED regional pseudobulk, at least three genes and 80% coverage. Scores and endpoints must have identical row order. The frozen seed is 20260829; patient bootstrap uses seed +101 and 1000 draws. All patient regions are held out together. CLI: python analysis/02_patient_external_transfer.py INPUT.npz --output transfer.json.
* 03 input is an NPZ with scores (122 x 201), labels (LE/IT/CT/PAN/MVP), groups (10 tumours). Frozen source construction averages present eligible member-gene z values; do not use rounded plotted scores. CLI: python analysis/03_anatomical_resolution.py INPUT.npz --output ivy.json. Association permutations use seed 2026083101, continuing the same RNG across ALL5, CORE3, PAN_MVP. Classification is the original nearest-centroid implementation. Permutation classification null is not included.
* 04 input is an NPZ with CARE_primary/CARE_recurrent (56 paired rows), GLASS_primary/GLASS_recurrent (111), GSE174554_primary/GSE174554_recurrent (20). Paired matrices must have identical patient and representation ordering and frozen complete-feature eligibility. Source scope is ALL_AUTHOR_ANNOTATED for CARE and GSE174554, WHOLE_TISSUE for GLASS. CLI: python analysis/04_longitudinal_deformation.py INPUT.npz --output longitudinal.json. It reproduces the LOPO all-k curve; supported dimension selection, finite-cell null, composition decomposition and mode gates are not included.
* 06 accepts a TSV column paired_difference, one row per patient. The frozen CRC04 preprocessing computes unsigned 11-gene segment means, patient-region medians, then invasive-front minus core. CLI: Rscript analysis/06_crc04_exact_inference.R INPUT.tsv. No argument runs a synthetic five-pair fixture only. The exact inference kernels preserve inclusive tails and zero handling. Real patient differences are not redistributed.

Never commit these prepared files or generated individual-level intermediates. Keep outputs outside the repository unless reviewed as safe aggregates. GBmap registry entries identify each author subset of the exact Core GBmap asset. GLASS identifies difg_glass,source data version2022-05-31. These sources are public, but exact prepared-score construction is not certified as Tier2.

## CRC04 public processed-input reproduction

Download Nanostring_GeoMx.zip from https://zenodo.org/records/17671259. Expected archive SHA256: da92403e7fe4381d1cc8c7f884ea71b7f461a3e880183924ccd5cfe5496681c4. Extract outside this repository and locate Nanostring_GeoMx/exp2/r_object/geomx_exp2_target_Data_norm.RDS. Expected object SHA256:5c4eea69b73a0244f0b475265dae93db0ca0806fd06606419c5de4d0bdab968a.

```sh
Rscript analysis/06_crc04_exact_inference.R --source-rds /path/to/geomx_exp2_target_Data_norm.RDS
```

The path above is a placeholder,not an investigator path. The module reads the source log_q_combat assay without renormalization,uses the original T1/PanCK+/CORE-or-INV selection,checks281 source segments and18441 features,selects85 eligible segments,and runs the frozen11-gene mean / patient-region median / exact sign-flip chain. It prints only aggregate results. The S4-slot access wrapper reads the same stored data as the original Biobase accessors without importing the entire GeoMx analysis stack.

## CARE source-dependent sampling and composition

```sh
python analysis/07_sampling_composition.py EXTERNAL_FROZEN_PROJECT_ROOT --output EXTERNAL_PRIVATE_OUTPUT
```

Input root must contain the frozen V1 context-stress, V2 context-deepening and V3 deformation directories named in the script. Required V3 prepared inputs are CARE_SPLIT_HALF_MODEL_SCORES.parquet, CARE_MATCHED_CELL_DOWNSAMPLE_MODEL_SCORES.parquet and CARE_BROAD_COMPARTMENT_REQUIRED_GENE_COUNTS.parquet. V1/V2 provide paired model scores, frozen gene mappings, pseudobulk counts and memberships. No external path is hard-coded. Output must be outside the public repository and frozen authority directories because intermediate outputs contain patient-level records. Only the aggregate summary is suitable for public certification. Seeds,1000 bootstrap draws,50/25-cell thresholds and the non-additive energy estimand are unchanged.

GBmap source version is independently identified by its deposited Core asset and local frozen collection metadata. See the source authors' [data availability](https://pmc.ncbi.nlm.nih.gov/articles/PMC12526130/) and the exact asset URL in datasets.tsv.
