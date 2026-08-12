suppressPackageStartupMessages({
  library(data.table)
  library(metafor)
  library(digest)
  library(jsonlite)
})

DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset="data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset="results")
REPO_ROOT <- Sys.getenv("GBM_MES_REPO_ROOT", unset=".")
RUN <- file.path(RESULT_ROOT,"spatial")
OUT <- file.path(RUN,"meta_test")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

patient_path <- file.path(RUN, "04_spatial_integration/SPATIAL_PATIENT_LEVEL.tsv")
meta_path <- file.path(RUN, "04_spatial_integration/SPATIAL_RANDOM_EFFECTS_META.tsv")
prediction_path <- file.path(RUN, "04_spatial_integration/SPATIAL_PREDICTION_INTERVAL.tsv")
final_path <- file.path(RUN, "FINAL_SPATIAL_V1_2.tsv")
script_path <- file.path(RUN, "scripts/spatial_v12_run.R")
protocol_path <- file.path(RUN, "00_governance/GBM_MES_SPATIAL_IVY_PROTOCOL_V1_2_FROZEN.md")
protocol_yaml_path <- file.path(RUN, "00_governance/freeze/protocol_freeze.yaml")
session_path <- file.path(RUN, "04_spatial_integration/SPATIAL_SESSION_INFO.txt")

required <- c(patient_path, meta_path, prediction_path, final_path, script_path,
              protocol_path, protocol_yaml_path, session_path)
stopifnot(all(file.exists(required)))

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE)
pat <- fread(patient_path)
frozen_meta <- fread(meta_path)
frozen_prediction <- fread(prediction_path)
final_spatial <- fread(final_path)

formal_rows <- vector("list", nrow(frozen_meta))
for (i in seq_len(nrow(frozen_meta))) {
  key <- frozen_meta[i]
  d <- pat[
    signature == key$signature & mode == key$mode & ecology == key$ecology &
      metric == key$metric & block_scale == key$block_scale & is.finite(patient_effect)
  ]
  x <- d$patient_effect
  ok <- is.finite(x) & abs(x) < 1
  x <- x[ok]
  slab <- paste(d$study[ok], d$patient_id[ok], sep = ":")
  z <- atanh(pmax(pmin(x, .999999), -.999999))
  vi <- rep(1 / 100, length(z))
  stopifnot(length(z) >= 2L)
  fit <- rma.uni(yi = z, vi = vi, method = "REML", slab = slab)
  pr <- predict(fit)
  formal_rows[[i]] <- data.table(
    signature = key$signature,
    mode = key$mode,
    ecology = key$ecology,
    metric = key$metric,
    block_scale = key$block_scale,
    study_n = uniqueN(d$study),
    patient_n = length(z),
    rho = tanh(as.numeric(fit$b)),
    ci_low = tanh(as.numeric(fit$ci.lb)),
    ci_high = tanh(as.numeric(fit$ci.ub)),
    prediction_low = tanh(as.numeric(pr$pi.lb)),
    prediction_high = tanh(as.numeric(pr$pi.ub)),
    tau2 = as.numeric(fit$tau2),
    i2 = as.numeric(fit$I2),
    test_statistic = as.numeric(fit$zval),
    df = NA_real_,
    p_value = as.numeric(fit$pval),
    formal_significance_status = ifelse(as.numeric(fit$pval) < 0.05, "P_LT_0_05", "P_GE_0_05"),
    formal_test_method = "metafor::rma.uni; Fisher-z patient effects; vi=0.01; REML tau2; intercept Wald z test; two-sided standard normal; test='z' default; df not applicable"
  )
}
formal <- rbindlist(formal_rows)

# The frozen family key is study x statistic x mode. At the cross-study level,
# study is the fixed CROSS_STUDY stratum; block scale is part of the statistic,
# matching the original block-q implementation. Each family therefore contains
# the prespecified 5 signatures x 4 ecology endpoints = 20 tests.
formal[, `:=`(
  fdr_method = "Benjamini-Hochberg",
  family_id = paste("CROSS_STUDY", metric, block_scale, mode, sep = "|"),
  family_size = .N,
  fdr = p.adjust(p_value, method = "BH")
), by = .(metric, block_scale, mode)]
stopifnot(all(formal$family_size == 20L))

pred_key <- c("signature", "mode", "ecology", "metric", "block_scale")
check <- merge(
  formal,
  frozen_meta,
  by = pred_key,
  suffixes = c("_reproduced", "_frozen"),
  all.x = TRUE,
  sort = FALSE
)
check <- merge(
  check,
  frozen_prediction[, c(pred_key, "prediction_low", "prediction_high"), with = FALSE],
  by = pred_key,
  suffixes = c("", "_frozen"),
  all.x = TRUE,
  sort = FALSE
)
setnames(check, c("prediction_low", "prediction_high"),
         c("prediction_low_reproduced", "prediction_high_reproduced"))

check[, `:=`(
  abs_error_rho = abs(rho - random_effect_rho),
  abs_error_ci_low = abs(ci_low_reproduced - ci_low_frozen),
  abs_error_ci_high = abs(ci_high_reproduced - ci_high_frozen),
  abs_error_prediction_low = abs(prediction_low_reproduced - prediction_low_frozen),
  abs_error_prediction_high = abs(prediction_high_reproduced - prediction_high_frozen),
  abs_error_tau2 = abs(tau2 - tau2_fisher_z),
  abs_error_i2 = abs(i2 - i2_percent)
)]
error_cols <- c("abs_error_rho", "abs_error_ci_low", "abs_error_ci_high",
                "abs_error_prediction_low", "abs_error_prediction_high",
                "abs_error_tau2", "abs_error_i2")
check[, max_abs_error := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = error_cols]
check[, reproduction_status := fifelse(max_abs_error <= 1e-10, "PASS_LE_1E_10", "FAIL_GT_1E_10")]

spot_check <- merge(
  formal[metric == "spot_spearman" & block_scale == "spot"],
  final_spatial,
  by = c("signature", "mode", "ecology"),
  suffixes = c("_reproduced", "_final"),
  all.x = TRUE,
  sort = FALSE
)
spot_check[, `:=`(
  abs_error_final_rho = abs(rho - patient_random_effect_rho),
  abs_error_final_ci_low = abs(ci_low_reproduced - ci_low_final),
  abs_error_final_ci_high = abs(ci_high_reproduced - ci_high_final),
  abs_error_final_prediction_low = abs(prediction_low_reproduced - prediction_low_final),
  abs_error_final_prediction_high = abs(prediction_high_reproduced - prediction_high_final),
  abs_error_final_i2 = abs(i2 - i2_percent_descriptive)
)]
spot_error_cols <- c("abs_error_final_rho", "abs_error_final_ci_low", "abs_error_final_ci_high",
                     "abs_error_final_prediction_low", "abs_error_final_prediction_high",
                     "abs_error_final_i2")
spot_check[, max_abs_error_final := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = spot_error_cols]
spot_check[, reproduction_status_final := fifelse(max_abs_error_final <= 1e-10, "PASS_LE_1E_10", "FAIL_GT_1E_10")]

overall_max <- max(c(check$max_abs_error, spot_check$max_abs_error_final), na.rm = TRUE)
reproduction_pass <- is.finite(overall_max) && overall_max <= 1e-10 &&
  all(check$reproduction_status == "PASS_LE_1E_10") &&
  all(spot_check$reproduction_status_final == "PASS_LE_1E_10")

check_out <- check[, c(
  pred_key,
  "rho", "random_effect_rho", "abs_error_rho",
  "ci_low_reproduced", "ci_low_frozen", "abs_error_ci_low",
  "ci_high_reproduced", "ci_high_frozen", "abs_error_ci_high",
  "prediction_low_reproduced", "prediction_low_frozen", "abs_error_prediction_low",
  "prediction_high_reproduced", "prediction_high_frozen", "abs_error_prediction_high",
  "tau2", "tau2_fisher_z", "abs_error_tau2",
  "i2", "i2_percent", "abs_error_i2", "max_abs_error", "reproduction_status"
), with = FALSE]

fwrite(formal[, .(
  signature, mode, ecology, metric, block_scale, study_n, patient_n,
  rho, ci_low, ci_high, prediction_low, prediction_high, tau2, i2,
  test_statistic, df, p_value, formal_significance_status, formal_test_method
)], file.path(OUT, "SPATIAL_RANDOM_EFFECTS_META_FORMAL_TEST_COMPLETION.tsv"), sep = "\t", na = "NA")

fwrite(formal[, .(
  signature, mode, ecology, metric, block_scale, study_n, patient_n,
  rho, ci_low, ci_high, prediction_low, prediction_high, tau2, i2,
  test_statistic, df, p_value, formal_significance_status, formal_test_method,
  fdr, fdr_method, family_id, family_size
)], file.path(OUT, "SPATIAL_RANDOM_EFFECTS_META_FDR_COMPLETION.tsv"), sep = "\t", na = "NA")

fwrite(formal[, .(
  family_id, fdr_method, family_size,
  signature, mode, ecology, metric, block_scale,
  family_definition = "CROSS_STUDY fixed stratum x statistic (metric+block_scale) x mode; all 5 prespecified signatures x 4 prespecified ecology endpoints",
  reference_source = protocol_yaml_path,
  reference_coordinate = "multiplicity.family_definition lines 97-100; original implementation spatial_v12_run.R lines 285-286"
)], file.path(OUT, "R6_MULTIPLICITY_FAMILY.tsv"), sep = "\t")

fwrite(check_out, file.path(OUT, "FROZEN_ESTIMATE_REPRODUCTION_CHECK.tsv"), sep = "\t")
fwrite(spot_check[, c(
  "signature", "mode", "ecology", "rho", "patient_random_effect_rho",
  "abs_error_final_rho", "ci_low_reproduced", "ci_low_final", "abs_error_final_ci_low",
  "ci_high_reproduced", "ci_high_final", "abs_error_final_ci_high",
  "prediction_low_reproduced", "prediction_low_final", "abs_error_final_prediction_low",
  "prediction_high_reproduced", "prediction_high_final", "abs_error_final_prediction_high",
  "i2", "i2_percent_descriptive", "abs_error_final_i2", "max_abs_error_final",
  "reproduction_status_final"
), with = FALSE], file.path(OUT, "FINAL_SPATIAL_V1_2_REPRODUCTION_CHECK.tsv"), sep = "\t")

method_lines <- c(
  "# R6 Meta Method Check",
  "",
  "Status: METHOD_RECOVERED",
  "",
  "This is a formal-test completion using frozen V1.2 patient-level summaries. No spatial score, gene set, overlap-pruning rule, patient/study inclusion rule, block scale, endpoint, or effect estimate was changed.",
  "",
  "## Exact recovered statistical specification",
  "",
  "- Effect input: patient-level Spearman rho, with sections equal-weighted within patient.",
  "- Effect scale: Fisher z = atanh(rho), clipped only at +/-0.999999 as in the original function.",
  "- Working variance: vi = 1/100 = 0.01 for every patient effect; equal-patient working variance, not spot-derived sampling variance.",
  "- Random-effects estimator: metafor::rma.uni(method = REML).",
  "- Formal test: intercept Wald z test; test='z' default; two-sided standard normal distribution.",
  "- Degrees of freedom: not applicable (NA); no t distribution and no Hartung-Knapp modification.",
  "- Confidence interval: metafor rma.uni default 95% Wald interval on Fisher-z scale, back-transformed with tanh.",
  "- Prediction interval: metafor::predict.rma default 95% normal-theory prediction interval, back-transformed with tanh.",
  "- Heterogeneity: REML tau^2 on Fisher-z scale; metafor I^2.",
  "- Multiplicity: Benjamini-Hochberg within CROSS_STUDY x statistic(metric+block_scale) x mode families; family size 20 (5 frozen signatures x 4 frozen ecology endpoints). RAW and OVERLAP_PRUNED are separate families.",
  "",
  "## Reference coordinates",
  "",
  paste0("- Exact model code: `", script_path, "`, lines 406-415; source SHA256 `", sha256(script_path), "`."),
  paste0("- Exact model inputs: `", patient_path, "`; source SHA256 `", sha256(patient_path), "`."),
  paste0("- Frozen meta target: `", meta_path, "`; source SHA256 `", sha256(meta_path), "`."),
  paste0("- Frozen prediction target: `", prediction_path, "`; source SHA256 `", sha256(prediction_path), "`."),
  paste0("- Final spot-level target: `", final_path, "`; source SHA256 `", sha256(final_path), "`."),
  paste0("- Protocol: `", protocol_path, "`, lines 24-26 and 42; source SHA256 `", sha256(protocol_path), "`."),
  paste0("- Machine-readable protocol: `", protocol_yaml_path, "`, multiplicity lines 97-100 and synthesis lines 101-107; source SHA256 `", sha256(protocol_yaml_path), "`."),
  paste0("- Original environment: `", session_path, "`; R 4.4.1 and metafor 4.8-0; source SHA256 `", sha256(session_path), "`."),
  "",
  "## Frozen estimate reproduction gate",
  "",
  paste0("- Rows reproduced against SPATIAL_RANDOM_EFFECTS_META.tsv: ", nrow(check), "."),
  paste0("- Spot rows reproduced against FINAL_SPATIAL_V1_2.tsv: ", nrow(spot_check), "."),
  paste0("- Maximum absolute numerical error: ", format(overall_max, digits = 17, scientific = TRUE), "."),
  paste0("- Required tolerance: <= 1e-10."),
  paste0("- Gate: ", ifelse(reproduction_pass, "PASS", "FAIL"), "."),
  "",
  "## Claim boundary",
  "",
  "Formal P/FDR completion does not upgrade the R6 claim or redefine program identities. MES1 remains BROAD_MES_ONLY; MES2 remains PARTIAL_FINE_GRAIN_RECOVERY; MES-Hyp remains SPATIAL_PATHOLOGY_ANCHORED; MES-Ast remains HETEROGENEOUS_ACROSS_SCALE."
)
writeLines(method_lines, file.path(OUT, "R6_META_METHOD_CHECK.md"))

status <- list(
  final_status = if (reproduction_pass) "FORMAL_TEST_COMPLETION_PASS" else "FROZEN_ESTIMATE_REPRODUCTION_FAIL",
  original_method_recovered = TRUE,
  frozen_estimates_reproduced = reproduction_pass,
  max_absolute_error = overall_max,
  tolerance = 1e-10,
  test_method = "REML Fisher-z random effects; Wald z; two-sided standard normal",
  fdr_method = "Benjamini-Hochberg",
  family_size = 20L,
  family_count = uniqueN(formal$family_id),
  formal_rows = nrow(formal)
)
write_json(status, file.path(OUT, "R6_FORMAL_TEST_COMPLETION_STATUS.json"), pretty = TRUE, auto_unbox = TRUE)

manifest_files <- list.files(OUT, full.names = TRUE)
manifest_files <- manifest_files[file.info(manifest_files)$isdir %in% FALSE]
manifest_files <- manifest_files[basename(manifest_files) != "R6_FORMAL_TEST_COMPLETION_MANIFEST_SHA256.tsv"]
manifest <- data.table(
  sha256 = vapply(manifest_files, sha256, character(1)),
  size_bytes = file.info(manifest_files)$size,
  file = basename(manifest_files)
)
fwrite(manifest, file.path(OUT, "R6_FORMAL_TEST_COMPLETION_MANIFEST_SHA256.tsv"), sep = "\t")

if (!reproduction_pass) stop("Frozen effect reproduction gate failed; formal completion not eligible for manuscript reference.")

cat(toJSON(status, pretty = TRUE, auto_unbox = TRUE), "\n")
