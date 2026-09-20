#!/usr/bin/env Rscript

DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset = "data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset = "results")
run_dir <- file.path(RESULT_ROOT, "care_realization")
care_source <- file.path(DATA_ROOT, "single_cell", "care", "source")
dir.create(file.path(run_dir, "work"), recursive = TRUE, showWarnings = FALSE)
axes <- c("MALIGNANT_STATE", "MYELOID", "HYPOXIA", "MATRIX", "VASCULAR_STROMAL")
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

multivar_fit <- function(y, x) {
  x <- as.matrix(x)
  design <- cbind(Intercept = 1, x)
  fit <- lm.fit(design, y)
  resid <- fit$residuals
  sse <- sum(resid^2)
  yc <- scale(y, center = TRUE, scale = FALSE)
  sst <- sum(yc^2)
  p <- ncol(x); n <- nrow(y)
  r2 <- 1 - sse / sst
  adj <- if (n - p - 1 > 0) 1 - (1 - r2) * (n - 1) / (n - p - 1) else NA_real_
  f <- if (p > 0) ((sst - sse) / p) / (sse / (n - p - 1)) else NA_real_
  list(sse = sse, sst = sst, r2 = r2, adjusted_r2 = adj, F = f, fitted = fit$fitted.values,
       residuals = resid, n = n, p = p)
}

omnibus_perm <- function(y, x, seed, b = 5000) {
  obs <- multivar_fit(y, x)
  set.seed(seed)
  null <- numeric(b)
  for (i in seq_len(b)) null[i] <- multivar_fit(y[sample(seq_len(nrow(y))), , drop = FALSE], x)$F
  c(F = obs$F, p = (1 + sum(null >= obs$F)) / (b + 1), null_median = median(null), null_P95 = quantile(null, .95))
}

partial_perm <- function(y, x_reduced, x_added, seed, b = 5000) {
  red <- multivar_fit(y, x_reduced)
  full <- multivar_fit(y, cbind(x_reduced, x_added))
  df_add <- ncol(as.matrix(x_added)); df_res <- nrow(y) - ncol(cbind(x_reduced, x_added)) - 1
  obs_f <- ((red$sse - full$sse) / df_add) / (full$sse / df_res)
  set.seed(seed)
  null <- numeric(b)
  for (i in seq_len(b)) {
    yp <- red$fitted + red$residuals[sample(seq_len(nrow(y))), , drop = FALSE]
    r <- multivar_fit(yp, x_reduced); f <- multivar_fit(yp, cbind(x_reduced, x_added))
    null[i] <- ((r$sse - f$sse) / df_add) / (f$sse / df_res)
  }
  c(F = obs_f, p = (1 + sum(null >= obs_f)) / (b + 1), null_median = median(null), null_P95 = quantile(null, .95))
}

driver <- as.data.frame(readRDS(file.path(care_source, "driver_snv_status_care_20250220.RDS")), stringsAsFactors = FALSE)
cna <- as.data.frame(readRDS(file.path(care_source, "genetic_select_gene_level_cna.RDS")), stringsAsFactors = FALSE)
arm <- as.data.frame(readRDS(file.path(care_source, "genetic_sample_level_chr_arm_cna.RDS")), stringsAsFactors = FALSE)
cell_meta <- as.data.frame(readRDS(file.path(care_source, "celltype_meta_data_2025_01_08.RDS")), stringsAsFactors = FALSE)
residual <- read.delim(file.path(run_dir, "P1_CARE_RESIDUAL_REALIZATION_MATRIX.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
scores <- read.delim(file.path(run_dir, "P1_CARE_PATIENT_LEVEL_AXES.tsv"), check.names = FALSE, stringsAsFactors = FALSE)

encode_genetics <- function(sample_ids) {
  out <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
  dmatch <- match(sample_ids, driver$ID)
  out$NF1_SNV <- driver$NF1[dmatch]
  out$EGFR_SNV <- driver$EGFR[dmatch]
  out$PTEN_SNV <- driver$PTEN[dmatch]
  cna_call <- function(gene, call) {
    vapply(sample_ids, function(s) {
      z <- cna[cna$ID == s & cna$gene_symbol == gene, "hlvl_call"]
      if (!length(z)) return(NA_real_)
      as.numeric(any(z == call, na.rm = TRUE))
    }, numeric(1))
  }
  out$EGFR_AMPLIFICATION <- cna_call("EGFR", 2)
  out$PDGFRA_AMPLIFICATION <- cna_call("PDGFRA", 2)
  out$PTEN_DEEP_DELETION <- cna_call("PTEN", -2)
  out$CDKN2A_DEEP_DELETION <- cna_call("CDKN2A", -2)
  arm_feature <- function(chr, desired) {
    vapply(sample_ids, function(s) {
      z <- arm[arm$ID == s & arm$chrom == chr, ]
      expected <- paste0(chr, c("p", "q"))
      vals <- setNames(z$arm_call, z$arm)[expected]
      if (any(is.na(vals))) return(NA_real_)
      as.numeric(all(vals == desired))
    }, numeric(1))
  }
  out$CHR7_GAIN <- arm_feature(7, 1)
  out$CHR10_LOSS <- arm_feature(10, -1)
  out
}

g_all <- encode_genetics(residual$sample_id)
candidate_names <- setdiff(names(g_all), "sample_id")
prevalence <- colSums(g_all[candidate_names] == 1, na.rm = TRUE)
available_n <- colSums(!is.na(g_all[candidate_names]))
feature_registry <- data.frame(
  feature = c("NF1_SNV", "EGFR_SNV", "EGFR_AMPLIFICATION", "PDGFRA_SNV", "PDGFRA_AMPLIFICATION",
              "PTEN_SNV", "PTEN_DEEP_DELETION", "CDKN2A_DEEP_DELETION", "CDKN2B_DEEP_DELETION",
              "CHR7_GAIN", "CHR10_LOSS", "ecDNA"),
  source = c(rep("driver_snv_status_care_20250220.RDS", 2), "genetic_select_gene_level_cna.RDS",
             "NOT_AVAILABLE_IN_OFFICIAL_SAMPLE_LEVEL_TABLE", "genetic_select_gene_level_cna.RDS",
             "driver_snv_status_care_20250220.RDS", rep("genetic_select_gene_level_cna.RDS", 3),
             rep("genetic_sample_level_chr_arm_cna.RDS", 2), "NO_FORMAL_PATIENT_LEVEL_ANNOTATION"),
  outcome_independent_selection = TRUE,
  stringsAsFactors = FALSE)
feature_registry$available_n <- vapply(feature_registry$feature, function(f) if (f %in% candidate_names) available_n[f] else 0, numeric(1))
feature_registry$altered_n <- vapply(feature_registry$feature, function(f) if (f %in% candidate_names) prevalence[f] else 0, numeric(1))
feature_registry$formal_eligibility <- feature_registry$available_n > 0 & feature_registry$altered_n >= 5
feature_registry$status <- ifelse(feature_registry$formal_eligibility, "FORMAL_PREVALENCE_GATE_PASS",
                                  ifelse(feature_registry$available_n > 0, "DESCRIPTIVE_PREVALENCE_LT5", "NOT_EVALUABLE_NO_FORMAL_CALL"))
feature_registry$definition <- c("SNV present", "SNV present", "high-level call +2", "not available", "high-level call +2",
                                 "SNV present", "deep deletion call -2", "deep deletion call -2; no CDKN2B sample-level row",
                                 "not available separately", "both 7p and 7q arm_call +1", "both 10p and 10q arm_call -1",
                                 "not inferred from CNA proxy")
formal_features <- intersect(feature_registry$feature[feature_registry$formal_eligibility], candidate_names)
# Reapply the frozen >=5-per-state rule within the exact multivariable
# complete-case set. This is outcome-independent and can only remove features.
repeat {
  complete_genetic <- complete.cases(g_all[, formal_features, drop = FALSE])
  exact <- g_all[complete_genetic, formal_features, drop = FALSE]
  exact_alt <- colSums(exact == 1)
  exact_ref <- colSums(exact == 0)
  keep <- exact_alt >= 5 & exact_ref >= 5
  if (all(keep)) break
  formal_features <- formal_features[keep]
  if (!length(formal_features)) stop("P3_NOT_EVALUABLE_NO_DRIVER_MEETS_EXACT_COMPLETE_CASE_PREVALENCE")
}
feature_registry$exact_joint_set_available_n <- ifelse(feature_registry$feature %in% formal_features, sum(complete_genetic), NA)
feature_registry$exact_joint_set_altered_n <- vapply(feature_registry$feature, function(f) {
  if (!f %in% formal_features) return(NA_real_)
  sum(g_all[complete_genetic, f] == 1)
}, numeric(1))
feature_registry$exact_joint_set_reference_n <- vapply(feature_registry$feature, function(f) {
  if (!f %in% formal_features) return(NA_real_)
  sum(g_all[complete_genetic, f] == 0)
}, numeric(1))
feature_registry$formal_multivariable_model <- feature_registry$feature %in% formal_features
feature_registry$status[feature_registry$formal_eligibility & !feature_registry$formal_multivariable_model] <- "DESCRIPTIVE_EXACT_JOINT_SET_STATE_LT5"
write_tsv(feature_registry, file.path(run_dir, "P3_GENOMIC_FEATURE_REGISTRY.tsv"))

# Non-tautological ecology block: compartments absent from the five-axis outcome definitions.
ct <- as.data.frame.matrix(table(cell_meta$ID, cell_meta$CellType))
for (nm in c("Malignant", "Lymphocyte", "Astrocyte", "Oligodendrocyte", "OPC", "Excitatory neuron", "Inhibitory neuron")) if (!nm %in% names(ct)) ct[[nm]] <- 0
ct$sample_id <- rownames(ct)
ct$total <- rowSums(ct[, setdiff(names(ct), "sample_id"), drop = FALSE])
ct$nonmalignant <- ct$total - ct$Malignant
ct$glioneural <- ct$Astrocyte + ct$Oligodendrocyte + ct$OPC + ct$`Excitatory neuron` + ct$`Inhibitory neuron`
ecology <- data.frame(sample_id = ct$sample_id,
                      lymphocyte_fraction_nonmalignant = ct$Lymphocyte / ct$nonmalignant,
                      glioneural_fraction_nonmalignant = ct$glioneural / ct$nonmalignant,
                      stringsAsFactors = FALSE)
ecology <- ecology[match(residual$sample_id, ecology$sample_id), ]
complete_ecology <- complete.cases(ecology[, -1, drop = FALSE])
analysis_keep <- complete_genetic & complete_ecology

y <- as.matrix(residual[analysis_keep, paste0("residual_z_", axes)])
g <- as.matrix(g_all[analysis_keep, formal_features, drop = FALSE])
e <- scale(as.matrix(ecology[analysis_keep, -1, drop = FALSE]))
rownames(y) <- residual$patient_id[analysis_keep]; rownames(g) <- rownames(y); rownames(e) <- rownames(y)
if (nrow(y) < 20 || ncol(g) < 1) stop("P3_NOT_EVALUABLE")

fit_g <- multivar_fit(y, g); fit_e <- multivar_fit(y, e); fit_ge <- multivar_fit(y, cbind(g, e))
perm_g <- omnibus_perm(y, g, 2026081831)
perm_e <- omnibus_perm(y, e, 2026081832)
perm_ge <- omnibus_perm(y, cbind(g, e), 2026081833)
partial_g <- partial_perm(y, e, g, 2026081833)
partial_e <- partial_perm(y, g, e, 2026081833)

fractions <- c(
  unique_genetics = fit_ge$adjusted_r2 - fit_e$adjusted_r2,
  unique_ecology = fit_ge$adjusted_r2 - fit_g$adjusted_r2,
  shared_genetics_ecology = fit_g$adjusted_r2 + fit_e$adjusted_r2 - fit_ge$adjusted_r2,
  unexplained = 1 - fit_ge$adjusted_r2
)

set.seed(2026081835)
boot_fraction <- matrix(NA_real_, 2000, 4, dimnames = list(NULL, names(fractions)))
for (rep in 1:2000) {
  ii <- sample(seq_len(nrow(y)), replace = TRUE)
  fg <- multivar_fit(y[ii, , drop = FALSE], g[ii, , drop = FALSE])
  fe <- multivar_fit(y[ii, , drop = FALSE], e[ii, , drop = FALSE])
  fge <- multivar_fit(y[ii, , drop = FALSE], cbind(g[ii, , drop = FALSE], e[ii, , drop = FALSE]))
  boot_fraction[rep, ] <- c(fge$adjusted_r2 - fe$adjusted_r2,
                            fge$adjusted_r2 - fg$adjusted_r2,
                            fg$adjusted_r2 + fe$adjusted_r2 - fge$adjusted_r2,
                            1 - fge$adjusted_r2)
}

partition <- data.frame(
  component = c("genetics_only_model", "ecology_only_model", "joint_model", names(fractions)),
  estimate = c(fit_g$adjusted_r2, fit_e$adjusted_r2, fit_ge$adjusted_r2, fractions),
  raw_R2 = c(fit_g$r2, fit_e$r2, fit_ge$r2, rep(NA, 4)),
  permutation_p = c(perm_g["p"], perm_e["p"], perm_ge["p"], partial_g["p"], partial_e["p"], NA, NA),
  test_role = c("OMNIBUS_GENETICS", "OMNIBUS_ECOLOGY", "OMNIBUS_JOINT", "PARTIAL_GENETICS_GIVEN_ECOLOGY",
                "PARTIAL_ECOLOGY_GIVEN_GENETICS", "ALGEBRAIC_SHARED_FRACTION", "ALGEBRAIC_UNEXPLAINED"),
  ci_low = c(rep(NA, 3), apply(boot_fraction, 2, quantile, .025)),
  ci_high = c(rep(NA, 3), apply(boot_fraction, 2, quantile, .975)),
  n_patients = nrow(y), genetic_features = paste(formal_features, collapse = ";"),
  ecology_features = "lymphocyte_fraction_nonmalignant;glioneural_fraction_nonmalignant",
  stringsAsFactors = FALSE)
write_tsv(partition, file.path(run_dir, "P3_VARIANCE_PARTITION.tsv"))

# Marginal rotation-invariant multivariate driver tests.
assoc <- list()
set.seed(2026081834)
for (j in seq_along(formal_features)) {
  x <- g[, j, drop = FALSE]
  obs <- multivar_fit(y, x)
  null <- numeric(5000)
  for (rep in 1:5000) null[rep] <- multivar_fit(y[sample(seq_len(nrow(y))), , drop = FALSE], x)$F
  p <- (1 + sum(null >= obs$F)) / 5001
  centroids <- aggregate(y, list(status = x[, 1]), mean)
  assoc[[j]] <- data.frame(feature = formal_features[j], n = nrow(y), altered_n = sum(x[, 1] == 1),
                           marginal_raw_R2 = obs$r2, marginal_adjusted_R2 = obs$adjusted_r2,
                           pseudo_F = obs$F, permutation_p = p,
                           centroid_distance = sqrt(sum((centroids[centroids$status == 1, -1] - centroids[centroids$status == 0, -1])^2)))
}
assoc <- do.call(rbind, assoc); assoc$BH_FDR <- p.adjust(assoc$permutation_p, "BH")
write_tsv(assoc, file.path(run_dir, "P3_REALIZATION_GENETIC_ASSOCIATIONS.tsv"))

# Longitudinal genomic stability sensitivity.
pair_registry <- read.delim(file.path(run_dir, "P2_PAIRED_PATIENT_REGISTRY.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
p2metrics <- read.delim(file.path(run_dir, "P2_REALIZATION_REWIRING_METRICS.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
ep <- pair_registry[pair_registry$eligibility, ]
gp <- encode_genetics(ep$primary_sample); gr <- encode_genetics(ep$recurrent_sample)
long_gen <- data.frame(patient_id = ep$patient_id, primary_sample = ep$primary_sample, recurrent_sample = ep$recurrent_sample,
                       genomics_complete = complete.cases(gp[, formal_features, drop = FALSE]) & complete.cases(gr[, formal_features, drop = FALSE]),
                       stringsAsFactors = FALSE)
long_gen$genomic_change_count <- NA_integer_; long_gen$genomic_stability_group <- "NOT_EVALUABLE"
for (i in seq_len(nrow(long_gen))) if (long_gen$genomics_complete[i]) {
  changes <- sum(as.numeric(gp[i, formal_features]) != as.numeric(gr[i, formal_features]))
  long_gen$genomic_change_count[i] <- changes
  long_gen$genomic_stability_group[i] <- if (changes == 0) "GENOMICALLY_STABLE" else "GENOMICALLY_ALTERED"
}
long_gen$realization_distance <- p2metrics$realization_distance_euclidean[match(long_gen$patient_id, p2metrics$patient_id)]
n_stable <- sum(long_gen$genomic_stability_group == "GENOMICALLY_STABLE")
n_altered <- sum(long_gen$genomic_stability_group == "GENOMICALLY_ALTERED")
formal_long <- n_stable >= 5 && n_altered >= 5
long_gen$comparison_status <- if (formal_long) "FORMAL_MIN5_PER_GROUP" else "DESCRIPTIVE_GROUP_TOO_SMALL"
if (formal_long) {
  wt <- wilcox.test(realization_distance ~ genomic_stability_group, data = long_gen[long_gen$genomic_stability_group != "NOT_EVALUABLE", ], exact = FALSE)
  long_gen$group_comparison_p <- wt$p.value
} else long_gen$group_comparison_p <- NA_real_
write_tsv(long_gen, file.path(run_dir, "P3_LONGITUDINAL_GENOMIC_SENSITIVITY.tsv"))

pg <- partial_g["p"]; pe <- partial_e["p"]
if (pg < .05 && pe < .05) {
  adjudication <- "GENETICS_AND_ECOLOGY_JOINTLY_CONSTRAIN_REALIZATION"
} else if (pe < .05 && fit_e$adjusted_r2 >= fit_g$adjusted_r2 + .05) {
  adjudication <- "ECOLOGICAL_CONTEXT_DOMINATES_REALIZATION_VARIANCE"
} else if (pg < .05 && fit_g$adjusted_r2 >= fit_e$adjusted_r2 + .05) {
  adjudication <- "GENETIC_ARCHITECTURE_DOMINATES_REALIZATION_VARIANCE"
} else if (perm_g["p"] < .05) {
  adjudication <- "GENETICS_CONSTRAINS_BUT_DOES_NOT_IDENTIFY_REALIZATION"
} else {
  adjudication <- "NO_STABLE_DETERMINANTS_IDENTIFIED"
}

md <- c("# P3 Determinant adjudication", "", paste0("Adjudication: `", adjudication, "`"), "",
        paste0("Complete patient-level determinant set: n=", nrow(y), "; formal genomic features: ", paste(formal_features, collapse = ", "), "."),
        paste0("Adjusted R2 genetics-only=", format(fit_g$adjusted_r2, digits = 4), " (permutation P=", format(perm_g["p"], digits = 4),
               "); ecology-only=", format(fit_e$adjusted_r2, digits = 4), " (P=", format(perm_e["p"], digits = 4),
               "); joint=", format(fit_ge$adjusted_r2, digits = 4), " (P=", format(perm_ge["p"], digits = 4), ")."),
        paste0("Variation partition: unique genetics=", format(fractions["unique_genetics"], digits = 4),
               ", unique ecology=", format(fractions["unique_ecology"], digits = 4),
               ", shared=", format(fractions["shared_genetics_ecology"], digits = 4),
               ", unexplained=", format(fractions["unexplained"], digits = 4), "."),
        "",
        "Ecology predictors deliberately exclude TAM, vascular/stromal, hypoxia and matrix variables used to define the realization outcome. The ecology fraction is therefore a restricted non-tautological estimate, not total ecological determination. Negative adjusted fractions are retained rather than truncated. Associations are constraints/explained variation, not causal effects.")
writeLines(md, file.path(run_dir, "P3_DETERMINANT_ADJUDICATION.md"))
saveRDS(list(y = y, genetics = g, ecology = e, partition = partition, associations = assoc,
             bootstrap_fraction = boot_fraction, adjudication = adjudication),
        file.path(run_dir, "work/P3_INTERNAL_RESULTS.rds"), compress = "xz")
message("P3_COMPLETE ", adjudication)
print(feature_registry); print(partition); print(assoc)
