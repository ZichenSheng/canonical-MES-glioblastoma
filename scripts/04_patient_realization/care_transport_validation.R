#!/usr/bin/env Rscript

DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset = "data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset = "results")
run_dir <- file.path(RESULT_ROOT, "care_realization")
care_input <- file.path(DATA_ROOT, "single_cell", "care", "prepared")
cgga_fingerprint <- file.path(RESULT_ROOT, "patient_realization", "02_fingerprint", "MES_CONDITIONED_FINGERPRINT.tsv")
dir.create(file.path(run_dir, "02_p1"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_dir, "work"), recursive = TRUE, showWarnings = FALSE)
axes <- c("MALIGNANT_STATE", "MYELOID", "HYPOXIA", "MATRIX", "VASCULAR_STROMAL")

write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
bh <- function(p) p.adjust(p, method = "BH")

pca_geometry <- function(x) {
  xc <- scale(as.matrix(x), center = TRUE, scale = FALSE)
  cv <- crossprod(xc) / (nrow(xc) - 1)
  eg <- eigen(cv, symmetric = TRUE)
  list(values = eg$values, vectors = eg$vectors, covariance = cv, scores = xc %*% eg$vectors)
}

subspace_metrics <- function(u, v) {
  s <- svd(t(u) %*% v, nu = 0, nv = 0)$d
  s <- sort(pmin(1, pmax(0, s)), decreasing = TRUE)
  angles <- acos(s) * 180 / pi
  k <- length(s)
  c(setNames(as.numeric(angles), paste0("theta", seq_len(k), "_degrees")),
    setNames(as.numeric(s), paste0("cos_theta", seq_len(k))),
    S_proj = sum(s^2) / k,
    D_chordal = sqrt(max(0, k - sum(s^2))))
}

score_from_reference <- function(expr, genes, baseline_samples) {
  use <- intersect(genes, rownames(expr))
  mu <- rowMeans(expr[use, baseline_samples, drop = FALSE])
  sig <- apply(expr[use, baseline_samples, drop = FALSE], 1, sd)
  keep <- is.finite(sig) & sig > 0
  use <- use[keep]; mu <- mu[keep]; sig <- sig[keep]
  z <- sweep(sweep(expr[use, , drop = FALSE], 1, mu, "-"), 1, sig, "/")
  list(score = colMeans(z), genes = use, mean = mu, sd = sig)
}

standardize_reference <- function(x, baseline_samples) {
  mu <- mean(x[baseline_samples]); sig <- sd(x[baseline_samples])
  if (!is.finite(sig) || sig <= 0) stop("Zero-variance patient-level axis")
  list(z = (x - mu) / sig, mean = mu, sd = sig)
}

pb <- readRDS(file.path(care_input, "CARE_TARGET_GENE_PSEUDOBULK_COUNTS.rds"))
registry <- read.delim(file.path(care_input, "CARE_SAMPLE_PATIENT_REGISTRY.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
elig <- read.delim(file.path(care_input, "CARE_ELIGIBILITY_AUDIT.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
registry <- merge(registry, elig[, c("sample_id", "sample_eligible", "p1_eligible")], by = "sample_id", all.x = TRUE, sort = FALSE)
registry <- registry[match(colnames(pb$all_counts), registry$sample_id), ]
stopifnot(identical(registry$sample_id, colnames(pb$all_counts)), sum(registry$p1_eligible) == 52L)
baseline_samples <- registry$sample_id[registry$p1_eligible]

all_cpm_log <- log2(sweep(pb$all_counts, 2, pb$all_library[colnames(pb$all_counts)] / 1e6, "/") + 1)
mal_cpm_log <- log2(sweep(pb$malignant_counts, 2, pb$malignant_library[colnames(pb$malignant_counts)] / 1e6, "/") + 1)

canonical_bulk_ref <- score_from_reference(all_cpm_log, pb$canonical, baseline_samples)
canonical_mal_ref <- score_from_reference(mal_cpm_log, pb$canonical, baseline_samples)
hypoxia_bulk_ref <- score_from_reference(all_cpm_log, pb$hypoxia_pruned, baseline_samples)
matrix_bulk_ref <- score_from_reference(all_cpm_log, pb$matrix_pruned, baseline_samples)
hypoxia_mal_ref <- score_from_reference(mal_cpm_log, pb$hypoxia_pruned, baseline_samples)
matrix_mal_ref <- score_from_reference(mal_cpm_log, pb$matrix_pruned, baseline_samples)

names(canonical_bulk_ref$score) <- colnames(all_cpm_log)
names(canonical_mal_ref$score) <- colnames(all_cpm_log)
names(hypoxia_bulk_ref$score) <- colnames(all_cpm_log)
names(matrix_bulk_ref$score) <- colnames(all_cpm_log)
names(hypoxia_mal_ref$score) <- colnames(all_cpm_log)
names(matrix_mal_ref$score) <- colnames(all_cpm_log)

myeloid_fraction <- setNames(registry$myeloid_cells / registry$total_evaluable_cells, registry$sample_id)
vascular_fraction <- setNames(registry$vascular_stromal_cells / registry$total_evaluable_cells, registry$sample_id)

b_std <- standardize_reference(canonical_bulk_ref$score, baseline_samples)
m_std <- standardize_reference(canonical_mal_ref$score, baseline_samples)
y_std <- standardize_reference(myeloid_fraction, baseline_samples)
h_std <- standardize_reference(hypoxia_bulk_ref$score, baseline_samples)
x_std <- standardize_reference(matrix_bulk_ref$score, baseline_samples)
v_std <- standardize_reference(vascular_fraction, baseline_samples)

score_table <- data.frame(
  patient_id = registry$patient_id,
  sample_id = registry$sample_id,
  timepoint = registry$timepoint,
  stage = registry$stage,
  sample_eligible = registry$sample_eligible,
  p1_cross_sectional_eligible = registry$p1_eligible,
  canonical_mes_raw = unname(canonical_bulk_ref$score[registry$sample_id]),
  canonical_mes_z = unname(b_std$z[registry$sample_id]),
  malignant_mes_raw = unname(canonical_mal_ref$score[registry$sample_id]),
  malignant_mes_z = unname(m_std$z[registry$sample_id]),
  myeloid_abundance_raw = unname(myeloid_fraction[registry$sample_id]),
  myeloid_z = unname(y_std$z[registry$sample_id]),
  hypoxia_raw = unname(hypoxia_bulk_ref$score[registry$sample_id]),
  hypoxia_z = unname(h_std$z[registry$sample_id]),
  matrix_raw = unname(matrix_bulk_ref$score[registry$sample_id]),
  matrix_z = unname(x_std$z[registry$sample_id]),
  vascular_stromal_abundance_raw = unname(vascular_fraction[registry$sample_id]),
  vascular_stromal_z = unname(v_std$z[registry$sample_id]),
  malignant_hypoxia_sensitivity_raw = unname(hypoxia_mal_ref$score[registry$sample_id]),
  malignant_matrix_sensitivity_raw = unname(matrix_mal_ref$score[registry$sample_id]),
  canonical_genes_effective = length(canonical_bulk_ref$genes),
  malignant_canonical_genes_effective = length(canonical_mal_ref$genes),
  hypoxia_pruned_genes_effective = length(hypoxia_bulk_ref$genes),
  matrix_pruned_genes_effective = length(matrix_bulk_ref$genes),
  stringsAsFactors = FALSE
)
write_tsv(score_table, file.path(run_dir, "P1_CARE_PATIENT_LEVEL_AXES.tsv"))

p1 <- score_table[score_table$p1_cross_sectional_eligible, ]
p1 <- p1[order(as.integer(sub("P", "", p1$patient_id))), ]
axis_cols <- c("malignant_mes_z", "myeloid_z", "hypoxia_z", "matrix_z", "vascular_stromal_z")
fmat <- as.matrix(p1[, axis_cols])
b <- p1$canonical_mes_z
residual_raw <- matrix(NA_real_, nrow(fmat), ncol(fmat), dimnames = list(p1$sample_id, axes))
for (k in seq_along(axes)) {
  for (i in seq_len(nrow(fmat))) {
    fit <- lm(fmat[-i, k] ~ b[-i])
    residual_raw[i, k] <- fmat[i, k] - (coef(fit)[1] + coef(fit)[2] * b[i])
  }
}
residual_z <- scale(residual_raw, center = TRUE, scale = TRUE)
residual_table <- data.frame(patient_id = p1$patient_id, sample_id = p1$sample_id, canonical_mes_z = b,
                             setNames(as.data.frame(residual_raw), paste0("residual_raw_", axes)),
                             setNames(as.data.frame(residual_z), paste0("residual_z_", axes)), check.names = FALSE)
write_tsv(residual_table, file.path(run_dir, "P1_CARE_RESIDUAL_REALIZATION_MATRIX.tsv"))

# Preserve full-baseline scoring and regression references for P2 projection.
full_regression <- lapply(seq_along(axes), function(k) lm(fmat[, k] ~ b))
names(full_regression) <- axes
full_res <- sapply(seq_along(axes), function(k) residuals(full_regression[[k]]))
colnames(full_res) <- axes
reference <- list(
  baseline_samples = p1$sample_id,
  gene_score_reference = list(canonical_bulk = canonical_bulk_ref, canonical_malignant = canonical_mal_ref,
                              hypoxia_bulk = hypoxia_bulk_ref, matrix_bulk = matrix_bulk_ref,
                              hypoxia_malignant = hypoxia_mal_ref, matrix_malignant = matrix_mal_ref),
  patient_axis_standardization = list(canonical = b_std, malignant = m_std, myeloid = y_std,
                                      hypoxia = h_std, matrix = x_std, vascular = v_std),
  residual_regression = full_regression,
  residual_mean = colMeans(full_res),
  residual_sd = apply(full_res, 2, sd),
  axes = axes
)
saveRDS(reference, file.path(run_dir, "work/P1_FROZEN_PROJECTION_REFERENCE.rds"), compress = "xz")

# Residual dependence and parallel-rank null share the exact frozen permutation stream.
corr_obs <- cor(residual_z, method = "spearman")
upper <- upper.tri(corr_obs)
t_obs <- sum(corr_obs[upper]^2)
set.seed(2026081801)
n_perm <- 5000L
t_null <- numeric(n_perm)
eig_null <- matrix(NA_real_, n_perm, 5)
for (rep in seq_len(n_perm)) {
  xp <- sapply(seq_len(5), function(j) sample(residual_z[, j], replace = FALSE))
  cm <- cor(xp, method = "spearman")
  t_null[rep] <- sum(cm[upper]^2)
  eig_null[rep, ] <- pca_geometry(xp)$values
}
global_p <- (1 + sum(t_null >= t_obs)) / (n_perm + 1)
global <- data.frame(n_patients = nrow(residual_z), statistic = "SUM_SQUARED_OFFDIAGONAL_SPEARMAN",
                     T_observed = t_obs, null_median = median(t_null), null_P95 = quantile(t_null, .95),
                     finite_permutation_p = global_p, permutations = n_perm,
                     status = if (global_p < .05) "STRUCTURED_DEPENDENCE_SUPPORTED" else "GLOBAL_DEPENDENCE_NOT_SUPPORTED")
write_tsv(global, file.path(run_dir, "P1_CARE_GLOBAL_DEPENDENCE_TEST.tsv"))

pair_idx <- which(upper, arr.ind = TRUE)
corr_rows <- lapply(seq_len(nrow(pair_idx)), function(i) {
  a <- pair_idx[i, 1]; d <- pair_idx[i, 2]
  test <- suppressWarnings(cor.test(residual_z[, a], residual_z[, d], method = "spearman", exact = FALSE))
  data.frame(axis_a = axes[a], axis_b = axes[d], spearman_rho = unname(test$estimate), nominal_p = test$p.value)
})
corr_df <- do.call(rbind, corr_rows); corr_df$BH_FDR <- bh(corr_df$nominal_p)
write_tsv(corr_df, file.path(run_dir, "02_p1/P1_CARE_RESIDUAL_CORRELATIONS.tsv"))

care_pca <- pca_geometry(residual_z)
variance_fraction <- care_pca$values / sum(care_pca$values)
parallel_p95 <- apply(eig_null, 2, quantile, .95)
effective_rank <- min(3L, sum(care_pca$values[1:3] > parallel_p95[1:3]))
rank_rows <- data.frame(
  component = paste0("PC", 1:5), eigenvalue = care_pca$values,
  variance_fraction = variance_fraction, cumulative_variance_fraction = cumsum(variance_fraction),
  permutation_eigenvalue_P95 = parallel_p95,
  exceeds_parallel_P95 = care_pca$values > parallel_p95,
  effective_rank_capped_3 = effective_rank,
  loading_MALIGNANT_STATE = care_pca$vectors[1, ], loading_MYELOID = care_pca$vectors[2, ],
  loading_HYPOXIA = care_pca$vectors[3, ], loading_MATRIX = care_pca$vectors[4, ],
  loading_VASCULAR_STROMAL = care_pca$vectors[5, ]
)
write_tsv(rank_rows, file.path(run_dir, "P1_CARE_RANK_ANALYSIS.tsv"))

rank_boot_summary <- list()
for (k in 1:3) {
  set.seed(2026081801 + k)
  ref <- care_pca$vectors[, 1:k, drop = FALSE]
  bmet <- matrix(NA_real_, 2000, k + 2)
  for (rep in 1:2000) {
    idx <- sample(seq_len(nrow(residual_z)), replace = TRUE)
    uu <- pca_geometry(residual_z[idx, , drop = FALSE])$vectors[, 1:k, drop = FALSE]
    met <- subspace_metrics(ref, uu)
    bmet[rep, 1] <- met["S_proj"]
    bmet[rep, 2] <- met["D_chordal"]
    bmet[rep, 2 + seq_len(k)] <- met[paste0("theta", seq_len(k), "_degrees")]
  }
  rank_boot_summary[[k]] <- data.frame(rank = k, replicates = 2000,
    median_S_proj = median(bmet[, 1]), S_proj_ci_low = quantile(bmet[, 1], .025), S_proj_ci_high = quantile(bmet[, 1], .975),
    median_D_chordal = median(bmet[, 2]),
    median_largest_angle_degrees = median(bmet[, k + 2]), largest_angle_P95_degrees = quantile(bmet[, k + 2], .95))
}
rank_boot <- do.call(rbind, rank_boot_summary)
write_tsv(rank_boot, file.path(run_dir, "02_p1/P1_CARE_RANK_BOOTSTRAP_SUMMARY.tsv"))

# CGGA bases reconstructed from the retained MES-conditioned fingerprint.
cgga <- read.delim(cgga_fingerprint, check.names = FALSE, stringsAsFactors = FALSE)
residual_cols <- paste0("residual_z_", axes)
cgga_bases <- lapply(c("CGGA325", "CGGA693"), function(cohort) {
  z <- as.matrix(cgga[cgga$cohort == cohort, residual_cols, drop = FALSE])
  pca_geometry(z)$vectors[, 1:2, drop = FALSE]
})
names(cgga_bases) <- c("CGGA325", "CGGA693")
ucare <- care_pca$vectors[, 1:2, drop = FALSE]

transport_rows <- list(); angle_rows <- list(); null_rows <- list()
for (anchor in names(cgga_bases)) {
  ua <- cgga_bases[[anchor]]
  obs <- subspace_metrics(ucare, ua)
  for (j in 1:2) angle_rows[[length(angle_rows) + 1L]] <- data.frame(
    comparison = paste0("CARE_VS_", anchor), principal_angle = paste0("theta", j),
    degrees = obs[paste0("theta", j, "_degrees")], canonical_correlation = obs[paste0("cos_theta", j)])

  set.seed(if (anchor == "CGGA325") 2026081805 else 2026081806)
  boot <- matrix(NA_real_, 2000, 3)
  colnames(boot) <- c("S_proj", "theta1", "theta2")
  for (rep in 1:2000) {
    idx <- sample(seq_len(nrow(residual_z)), replace = TRUE)
    ub <- pca_geometry(residual_z[idx, , drop = FALSE])$vectors[, 1:2, drop = FALSE]
    bm <- subspace_metrics(ub, ua)
    boot[rep, ] <- bm[c("S_proj", "theta1_degrees", "theta2_degrees")]
  }

  perms <- as.matrix(expand.grid(1:5, 1:5, 1:5, 1:5, 1:5))
  perms <- perms[apply(perms, 1, function(z) length(unique(z)) == 5), , drop = FALSE]
  exact_s <- apply(perms, 1, function(pp) subspace_metrics(ucare[pp, , drop = FALSE], ua)["S_proj"])
  exact_p <- mean(exact_s >= obs["S_proj"] - 1e-15)

  set.seed(2026081807)
  ind_s <- numeric(5000)
  for (rep in 1:5000) {
    xp <- sapply(seq_len(5), function(j) sample(residual_z[, j], replace = FALSE))
    up <- pca_geometry(xp)$vectors[, 1:2, drop = FALSE]
    ind_s[rep] <- subspace_metrics(up, ua)["S_proj"]
  }
  ind_p <- (1 + sum(ind_s >= obs["S_proj"] - 1e-15)) / 5001

  set.seed(2026081808)
  grass_s <- numeric(100000)
  for (rep in 1:100000) {
    q <- qr.Q(qr(matrix(rnorm(10), 5, 2)))[, 1:2, drop = FALSE]
    grass_s[rep] <- subspace_metrics(q, ua)["S_proj"]
  }
  grass_p <- mean(grass_s >= obs["S_proj"] - 1e-15)

  transport_pass <- obs["S_proj"] >= .80 && obs["theta2_degrees"] <= 25 &&
    ind_p <= .01 && quantile(boot[, "S_proj"], .025) >= .65
  transport_rows[[anchor]] <- data.frame(
    comparison = paste0("CARE_VS_", anchor), n_CARE = nrow(residual_z),
    theta1_degrees = obs["theta1_degrees"], theta2_degrees = obs["theta2_degrees"],
    canonical_correlation_1 = obs["cos_theta1"], canonical_correlation_2 = obs["cos_theta2"],
    projection_similarity = obs["S_proj"], grassmann_chordal_distance = obs["D_chordal"],
    bootstrap_median_S_proj = median(boot[, "S_proj"]),
    bootstrap_S_proj_ci_low = quantile(boot[, "S_proj"], .025),
    bootstrap_S_proj_ci_high = quantile(boot[, "S_proj"], .975),
    exact_axis_label_p = exact_p, independent_axis_permutation_p = ind_p,
    grassmann_geometric_tail = grass_p, transport_gate_pass = transport_pass)
  null_rows[[length(null_rows) + 1L]] <- data.frame(comparison = paste0("CARE_VS_", anchor), null_model = "EXACT_AXIS_LABEL_PERMUTATION", iterations = 120, observed_S_proj = obs["S_proj"], null_median = median(exact_s), null_P95 = quantile(exact_s, .95), tail_probability = exact_p)
  null_rows[[length(null_rows) + 1L]] <- data.frame(comparison = paste0("CARE_VS_", anchor), null_model = "INDEPENDENT_AXIS_PATIENT_PERMUTATION", iterations = 5000, observed_S_proj = obs["S_proj"], null_median = median(ind_s), null_P95 = quantile(ind_s, .95), tail_probability = ind_p)
  null_rows[[length(null_rows) + 1L]] <- data.frame(comparison = paste0("CARE_VS_", anchor), null_model = "RANDOM_GRASSMANN_GEOMETRY", iterations = 100000, observed_S_proj = obs["S_proj"], null_median = median(grass_s), null_P95 = quantile(grass_s, .95), tail_probability = grass_p)
}

transport <- do.call(rbind, transport_rows)
angles <- do.call(rbind, angle_rows)
nulls <- do.call(rbind, null_rows)
write_tsv(transport, file.path(run_dir, "P1_CARE_RANK2_TRANSPORT.tsv"))
write_tsv(angles, file.path(run_dir, "P1_CARE_PRINCIPAL_ANGLES.tsv"))
write_tsv(nulls, file.path(run_dir, "P1_CARE_SUBSPACE_NULL.tsv"))

rank2_internal <- rank_boot[rank_boot$rank == 2, ]
internal_pass <- rank2_internal$median_S_proj >= .80 && rank2_internal$S_proj_ci_low >= .60 && rank2_internal$median_largest_angle_degrees <= 25
if (global_p >= .05) {
  adjudication <- "CELL_RESOLVED_REALIZATION_STRUCTURE_NOT_REPLICATED"
} else if (effective_rank != 2L) {
  adjudication <- "CELL_RESOLVED_STRUCTURED_REALIZATION_SUPPORTED_RANK_NOT_IDENTICAL"
} else if (internal_pass && all(transport$transport_gate_pass)) {
  adjudication <- "CELL_RESOLVED_RANK2_TRANSPORT_SUPPORTED"
} else {
  adjudication <- "STRUCTURE_PRESENT_BUT_SUBSPACE_NOT_TRANSPORTED"
}

md <- c(
  "# P1 Cell-resolved transport adjudication",
  "",
  paste0("Adjudication: `", adjudication, "`"),
  "",
  paste0("Eligible patient-level CARE ecosystems: ", nrow(residual_z), "."),
  paste0("Global residual-dependence finite-permutation P: ", format(global_p, digits = 5), "."),
  paste0("Prespecified CARE effective rank (parallel null; capped at 3): ", effective_rank, "."),
  paste0("CARE rank-2 cumulative variance: ", format(sum(variance_fraction[1:2]), digits = 5), "."),
  "",
  paste0("CARE–CGGA325 S_proj=", format(transport$projection_similarity[transport$comparison == "CARE_VS_CGGA325"], digits = 5),
         "; angles=", format(transport$theta1_degrees[transport$comparison == "CARE_VS_CGGA325"], digits = 4), "°, ",
         format(transport$theta2_degrees[transport$comparison == "CARE_VS_CGGA325"], digits = 4), "°."),
  paste0("CARE–CGGA693 S_proj=", format(transport$projection_similarity[transport$comparison == "CARE_VS_CGGA693"], digits = 5),
         "; angles=", format(transport$theta1_degrees[transport$comparison == "CARE_VS_CGGA693"], digits = 4), "°, ",
         format(transport$theta2_degrees[transport$comparison == "CARE_VS_CGGA693"], digits = 4), "°."),
  "",
  "This is a rotation-invariant frozen-construct transport test. It does not validate BayesPrism as an algorithm, establish universal rank two, define new MES axes, or support causal interpretation."
)
writeLines(md, file.path(run_dir, "P1_CELL_RESOLVED_TRANSPORT_ADJUDICATION.md"))

saveRDS(list(residual_z = residual_z, pca = care_pca, t_null = t_null, eigen_null = eig_null,
             rank_bootstrap = rank_boot, transport = transport, adjudication = adjudication),
        file.path(run_dir, "work/P1_INTERNAL_RESULTS.rds"), compress = "xz")
message("P1_COMPLETE ", adjudication)
print(global)
print(rank_rows[, 1:7])
print(transport)
