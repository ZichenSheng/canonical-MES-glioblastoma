#!/usr/bin/env Rscript

DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset = "data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset = "results")
run_dir <- file.path(RESULT_ROOT, "care_realization")
care_input <- file.path(DATA_ROOT, "single_cell", "care", "prepared")
dir.create(file.path(run_dir, "03_p2"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_dir, "work"), recursive = TRUE, showWarnings = FALSE)
axes <- c("MALIGNANT_STATE", "MYELOID", "HYPOXIA", "MATRIX", "VASCULAR_STROMAL")
axis_cols <- c("malignant_mes_z", "myeloid_z", "hypoxia_z", "matrix_z", "vascular_stromal_z")
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

scores <- read.delim(file.path(run_dir, "P1_CARE_PATIENT_LEVEL_AXES.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
registry <- read.delim(file.path(care_input, "CARE_SAMPLE_PATIENT_REGISTRY.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
ref <- readRDS(file.path(run_dir, "work/P1_FROZEN_PROJECTION_REFERENCE.rds"))
p1res <- read.delim(file.path(run_dir, "P1_CARE_RESIDUAL_REALIZATION_MATRIX.tsv"), check.names = FALSE, stringsAsFactors = FALSE)

fmat_all <- as.matrix(scores[, axis_cols])
b_all <- scores$canonical_mes_z
res_all <- matrix(NA_real_, nrow(scores), 5, dimnames = list(scores$sample_id, axes))
for (k in seq_along(axes)) {
  cf <- coef(ref$residual_regression[[axes[k]]])
  raw <- fmat_all[, k] - (cf[1] + cf[2] * b_all)
  res_all[, k] <- (raw - ref$residual_mean[k]) / ref$residual_sd[k]
}

time_num <- function(x) as.integer(sub("T", "", x))
patient_ids <- sort(unique(registry$patient_id), index.return = FALSE)
pair_rows <- list()
for (pid in patient_ids) {
  z <- registry[registry$patient_id == pid, ]
  prim <- z[z$stage == "Primary", ]
  recur <- z[grepl("Recurrent", z$stage), ]
  if (!nrow(prim) || !nrow(recur)) {
    pair_rows[[length(pair_rows) + 1L]] <- data.frame(patient_id = pid, primary_sample = NA, recurrent_sample = NA,
      interval_month = NA, primary_cells = NA, recurrent_cells = NA, primary_malignant_cells = NA,
      recurrent_malignant_cells = NA, treatment_history = NA, axis_completeness = FALSE,
      eligibility = FALSE, exclusion_reason = "NO_CERTIFIED_PRIMARY_TO_RECURRENCE_PAIR")
    next
  }
  prim <- prim[which.min(time_num(prim$timepoint)), , drop = FALSE]
  recur <- recur[time_num(recur$timepoint) > time_num(prim$timepoint), , drop = FALSE]
  if (!nrow(recur)) {
    pair_rows[[length(pair_rows) + 1L]] <- data.frame(patient_id = pid, primary_sample = prim$sample_id, recurrent_sample = NA,
      interval_month = NA, primary_cells = prim$total_evaluable_cells, recurrent_cells = NA,
      primary_malignant_cells = prim$malignant_cells, recurrent_malignant_cells = NA, treatment_history = NA,
      axis_completeness = FALSE, eligibility = FALSE, exclusion_reason = "RECURRENCE_NOT_LATER_THAN_PRIMARY")
    next
  }
  recur <- recur[which.min(time_num(recur$timepoint)), , drop = FALSE]
  ps <- scores[scores$sample_id == prim$sample_id, ]; rs <- scores[scores$sample_id == recur$sample_id, ]
  complete <- nrow(ps) == 1L && nrow(rs) == 1L && all(is.finite(c(ps$canonical_mes_z, unlist(ps[axis_cols]), rs$canonical_mes_z, unlist(rs[axis_cols]))))
  eligible <- isTRUE(ps$sample_eligible) && isTRUE(rs$sample_eligible) && complete
  reason <- if (eligible) "NONE" else paste0(if (!isTRUE(ps$sample_eligible)) "PRIMARY_SAMPLE_GATE;" else "",
                                              if (!isTRUE(rs$sample_eligible)) "RECURRENT_SAMPLE_GATE;" else "",
                                              if (!complete) "AXIS_INCOMPLETE;" else "")
  treatment <- paste0("radiation_before_recurrence=", recur$radiation_before_surgery,
                      ";alkylate_before_recurrence=", recur$alkylate_before_surgery,
                      ";steroid_before_recurrence=", recur$steroid_before_surgery)
  pair_rows[[length(pair_rows) + 1L]] <- data.frame(
    patient_id = pid, primary_sample = prim$sample_id, recurrent_sample = recur$sample_id,
    interval_month = recur$surgical_interval_month,
    primary_cells = prim$total_evaluable_cells, recurrent_cells = recur$total_evaluable_cells,
    primary_malignant_cells = prim$malignant_cells, recurrent_malignant_cells = recur$malignant_cells,
    treatment_history = treatment, axis_completeness = complete, eligibility = eligible, exclusion_reason = reason,
    stringsAsFactors = FALSE)
}
pair_registry <- do.call(rbind, pair_rows)
write_tsv(pair_registry, file.path(run_dir, "P2_PAIRED_PATIENT_REGISTRY.tsv"))

eligible_pairs <- pair_registry[pair_registry$eligibility, ]
if (nrow(eligible_pairs) < 15) stop("NOT_EVALUABLE_DUE_TO_PAIRING_OR_POWER")
stopifnot(nrow(eligible_pairs) == 42L, !anyDuplicated(eligible_pairs$patient_id))

long_rows <- list(); metric_rows <- list()
base_cov <- cov(as.matrix(p1res[, paste0("residual_z_", axes)]))
base_cov_inv <- solve(base_cov)
delta_matrix <- matrix(NA_real_, nrow(eligible_pairs), 5, dimnames = list(eligible_pairs$patient_id, axes))
primary_matrix <- recurrent_matrix <- delta_matrix
primary_mes <- recurrent_mes <- setNames(numeric(nrow(eligible_pairs)), eligible_pairs$patient_id)
for (i in seq_len(nrow(eligible_pairs))) {
  pr <- eligible_pairs$primary_sample[i]; rr <- eligible_pairs$recurrent_sample[i]; pid <- eligible_pairs$patient_id[i]
  ip <- match(pr, scores$sample_id); ir <- match(rr, scores$sample_id)
  rp <- res_all[pr, ]; rrvec <- res_all[rr, ]; d <- rrvec - rp
  primary_matrix[pid, ] <- rp; recurrent_matrix[pid, ] <- rrvec; delta_matrix[pid, ] <- d
  primary_mes[pid] <- scores$canonical_mes_z[ip]; recurrent_mes[pid] <- scores$canonical_mes_z[ir]
  for (tp in c("PRIMARY", "RECURRENT")) {
    sid <- if (tp == "PRIMARY") pr else rr
    ii <- if (tp == "PRIMARY") ip else ir
    rv <- if (tp == "PRIMARY") rp else rrvec
    long_rows[[length(long_rows) + 1L]] <- data.frame(patient_id = pid, timepoint_role = tp, sample_id = sid,
      stage = scores$stage[ii], canonical_mes_z = scores$canonical_mes_z[ii],
      malignant_mes_z = scores$malignant_mes_z[ii], myeloid_z = scores$myeloid_z[ii],
      hypoxia_z = scores$hypoxia_z[ii], matrix_z = scores$matrix_z[ii], vascular_stromal_z = scores$vascular_stromal_z[ii],
      residual_MALIGNANT_STATE = rv[1], residual_MYELOID = rv[2], residual_HYPOXIA = rv[3],
      residual_MATRIX = rv[4], residual_VASCULAR_STROMAL = rv[5])
  }
  euclidean <- sqrt(sum(d^2)); mahal <- sqrt(as.numeric(t(d) %*% base_cov_inv %*% d))
  metric_rows[[i]] <- data.frame(patient_id = pid, primary_sample = pr, recurrent_sample = rr,
    delta_MES = recurrent_mes[pid] - primary_mes[pid], abs_delta_MES = abs(recurrent_mes[pid] - primary_mes[pid]),
    delta_MALIGNANT_STATE = d[1], delta_MYELOID = d[2], delta_HYPOXIA = d[3], delta_MATRIX = d[4], delta_VASCULAR_STROMAL = d[5],
    realization_distance_euclidean = euclidean, realization_distance_mahalanobis_sensitivity = mahal)
}
longitudinal <- do.call(rbind, long_rows)
metrics <- do.call(rbind, metric_rows)
median_distance <- median(metrics$realization_distance_euclidean)
metrics$illustrative_iso_MES <- metrics$abs_delta_MES <= .25
metrics$illustrative_high_rewiring <- metrics$realization_distance_euclidean >= median_distance
metrics$illustrative_iso_MES_high_rewiring <- metrics$illustrative_iso_MES & metrics$illustrative_high_rewiring
write_tsv(longitudinal, file.path(run_dir, "P2_LONGITUDINAL_REALIZATION_MATRIX.tsv"))

rho_test <- suppressWarnings(cor.test(metrics$abs_delta_MES, metrics$realization_distance_euclidean, method = "spearman", exact = FALSE))
fit <- lm(realization_distance_euclidean ~ abs_delta_MES, data = metrics)
metrics$distance_unexplained_by_abs_delta_MES <- residuals(fit)
write_tsv(metrics, file.path(run_dir, "P2_REALIZATION_REWIRING_METRICS.tsv"))

# Pair-level bootstrap for primary continuous summaries.
set.seed(2026081821)
boot <- matrix(NA_real_, 2000, 4, dimnames = list(NULL, c("median_distance", "spearman_rho", "ols_R2", "iso_high_fraction")))
for (rep in 1:2000) {
  ii <- sample(seq_len(nrow(metrics)), replace = TRUE)
  z <- metrics[ii, ]
  boot[rep, "median_distance"] <- median(z$realization_distance_euclidean)
  boot[rep, "spearman_rho"] <- suppressWarnings(cor(z$abs_delta_MES, z$realization_distance_euclidean, method = "spearman"))
  boot[rep, "ols_R2"] <- summary(lm(realization_distance_euclidean ~ abs_delta_MES, data = z))$r.squared
  boot[rep, "iso_high_fraction"] <- mean(z$illustrative_iso_MES_high_rewiring)
}
summary_table <- data.frame(
  metric = c("median_realization_distance", "spearman_abs_deltaMES_vs_distance", "OLS_R2_distance_on_abs_deltaMES", "illustrative_iso_MES_high_rewiring_fraction"),
  estimate = c(median_distance, unname(rho_test$estimate), summary(fit)$r.squared, mean(metrics$illustrative_iso_MES_high_rewiring)),
  ci_low = apply(boot, 2, quantile, .025), ci_high = apply(boot, 2, quantile, .975),
  nominal_p = c(NA, rho_test$p.value, coef(summary(fit))[2, 4], NA))
write_tsv(summary_table, file.path(run_dir, "03_p2/P2_CONTINUOUS_SUMMARY.tsv"))

# Null 1: pair permutation, one-to-one recurrent reassignment.
set.seed(2026081822)
pair_null <- numeric(5000)
for (rep in 1:5000) {
  pp <- sample(seq_len(nrow(metrics)), replace = FALSE)
  dd <- recurrent_matrix[pp, , drop = FALSE] - primary_matrix
  pair_null[rep] <- median(sqrt(rowSums(dd^2)))
}
pair_p <- (1 + sum(pair_null >= median_distance)) / 5001

# Null 2: signed time-label shuffle; magnitude is algebraically invariant.
obs_direction <- sqrt(sum(colMeans(delta_matrix)^2))
set.seed(2026081823)
direction_null <- numeric(5000)
for (rep in 1:5000) {
  signs <- sample(c(-1, 1), nrow(delta_matrix), replace = TRUE)
  direction_null[rep] <- sqrt(sum(colMeans(delta_matrix * signs)^2))
}
direction_p <- (1 + sum(direction_null >= obs_direction)) / 5001

# Null 3: non-self recurrence matched to each observed absolute deltaMES.
set.seed(2026081824)
matched_null <- numeric(5000)
for (rep in 1:5000) {
  dists <- numeric(nrow(metrics))
  for (i in seq_len(nrow(metrics))) {
    candidates <- setdiff(seq_len(nrow(metrics)), i)
    candidate_amp <- abs(recurrent_mes[candidates] - primary_mes[i])
    mismatch <- abs(candidate_amp - metrics$abs_delta_MES[i])
    best <- candidates[mismatch == min(mismatch)]
    j <- sample(best, 1)
    dists[i] <- sqrt(sum((recurrent_matrix[j, ] - primary_matrix[i, ])^2))
  }
  matched_null[rep] <- median(dists)
}
matched_p <- (1 + sum(matched_null >= median_distance)) / 5001

null_table <- data.frame(
  null_model = c("PAIR_PERMUTATION", "TIMEPOINT_LABEL_SHUFFLE_DIRECTION_ONLY", "DELTA_MES_MATCHED_NONSELF_RECURRENCE"),
  statistic = c("median_euclidean_distance", "norm_mean_deltaR", "median_euclidean_distance"),
  observed = c(median_distance, obs_direction, median_distance),
  null_median = c(median(pair_null), median(direction_null), median(matched_null)),
  null_P95 = c(quantile(pair_null, .95), quantile(direction_null, .95), quantile(matched_null, .95)),
  upper_tail_finite_p = c(pair_p, direction_p, matched_p), iterations = 5000,
  interpretation = c("Random recurrent patient pairing", "Distance magnitude invariant; direction test only", "Random non-self recurrence matched to observed abs(deltaMES)"))
write_tsv(null_table, file.path(run_dir, "P2_ISO_MES_NULL_TESTS.tsv"))

# Signed axis changes are secondary and multiplicity-controlled.
axis_change <- lapply(seq_along(axes), function(k) {
  wt <- suppressWarnings(wilcox.test(delta_matrix[, k], mu = 0, exact = FALSE))
  data.frame(axis = axes[k], median_delta = median(delta_matrix[, k]), mean_delta = mean(delta_matrix[, k]),
             positive_fraction = mean(delta_matrix[, k] > 0), wilcoxon_p = wt$p.value)
})
axis_change <- do.call(rbind, axis_change); axis_change$BH_FDR <- p.adjust(axis_change$wilcoxon_p, "BH")
write_tsv(axis_change, file.path(run_dir, "03_p2/P2_SIGNED_AXIS_CHANGES.tsv"))

parallel_orthogonal <- data.frame(
  patient_id = metrics$patient_id, delta_MES = metrics$delta_MES,
  MES_parallel_component = NA_real_, MES_orthogonal_component = NA_real_,
  residual_realization_distance = metrics$realization_distance_euclidean,
  status = "NOT_EVALUABLE_NO_FROZEN_COMMON_MES_DIRECTION",
  reason = "The retained analysis contains axis-specific MES residual regressions but no prespecified common five-axis MES direction; no post hoc direction was introduced.")
write_tsv(parallel_orthogonal, file.path(run_dir, "P2_PARALLEL_ORTHOGONAL_CHANGE.tsv"))

iso_n <- sum(metrics$illustrative_iso_MES_high_rewiring)
support <- pair_p < .05 && matched_p < .05 && summary(fit)$r.squared < .25 && abs(unname(rho_test$estimate)) < .50 && iso_n >= 5 && iso_n / nrow(metrics) >= .10
if (support) {
  adjudication <- "LONGITUDINAL_ISO_MES_REWIRING_SUPPORTED"
} else if (summary(fit)$r.squared >= .25 || abs(unname(rho_test$estimate)) >= .50) {
  adjudication <- "LONGITUDINAL_REALIZATION_REWIRING_PRESENT_BUT_MES_DEPENDENT"
} else if (IQR(metrics$realization_distance_euclidean) > 0 && direction_p >= .05) {
  adjudication <- "HETEROGENEOUS_PATIENT_SPECIFIC_REWIRING"
} else {
  adjudication <- "NO_REPRODUCIBLE_REWIRING"
}

md <- c(
  "# P2 Longitudinal rewiring adjudication", "", paste0("Adjudication: `", adjudication, "`"), "",
  paste0("Strict eligible primary-to-earliest-recurrence pairs: ", nrow(metrics), "."),
  paste0("Median z-standardized Euclidean realization change: ", format(median_distance, digits = 5),
         " (95% pair-bootstrap CI ", format(quantile(boot[, "median_distance"], .025), digits = 4), "–",
         format(quantile(boot[, "median_distance"], .975), digits = 4), ")."),
  paste0("Spearman abs(deltaMES) versus realization distance: rho=", format(unname(rho_test$estimate), digits = 4),
         ", P=", format(rho_test$p.value, digits = 4), "; OLS R2=", format(summary(fit)$r.squared, digits = 4), "."),
  paste0("Pair-permutation upper-tail P=", format(pair_p, digits = 5), "; deltaMES-matched upper-tail P=", format(matched_p, digits = 5), "."),
  paste0("Illustrative abs(deltaMES)<=0.25 plus distance>=cohort median: ", iso_n, "/", nrow(metrics), " pairs. This subgroup is descriptive only."),
  "",
  "Time-label shuffling cannot change Euclidean distance; it was correctly used only for signed direction. Results describe disease-course-associated reconfiguration and do not identify a treatment effect."
)
writeLines(md, file.path(run_dir, "P2_LONGITUDINAL_REWIRING_ADJUDICATION.md"))
saveRDS(list(metrics = metrics, delta_matrix = delta_matrix, pair_null = pair_null, direction_null = direction_null,
             matched_null = matched_null, continuous_bootstrap = boot, adjudication = adjudication),
        file.path(run_dir, "work/P2_INTERNAL_RESULTS.rds"), compress = "xz")
message("P2_COMPLETE ", adjudication)
print(summary_table); print(null_table); print(axis_change)
