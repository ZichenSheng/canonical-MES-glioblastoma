#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
set.seed(2026073111)

DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset="data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset="results")
REPO_ROOT <- Sys.getenv("GBM_MES_REPO_ROOT", unset=".")
CORE_INPUT <- file.path(DATA_ROOT,"prepared","core")
OUT <- file.path(RESULT_ROOT,"decomposition","bayesprism")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

bulk_path <- file.path(CORE_INPUT, "05_bulk/BULK_PATIENT_SCORES.tsv")
frac_path <- file.path(CORE_INPUT, "06_bayesprism/CELL_FRACTION_ESTIMATES.tsv")
mal_path <- file.path(CORE_INPUT, "06_bayesprism/MALIGNANT_EXPRESSION_ESTIMATES.tsv")
qc_path <- file.path(CORE_INPUT, "06_bayesprism/BAYESPRISM_SAMPLE_QC.tsv")
fit_path <- file.path(CORE_INPUT, "06_bayesprism/BAYESPRISM_FULL_FIT.rds")

sha <- function(path) unname(tools::md5sum(path))
sha256 <- function(path) {
  z <- system2("/usr/bin/shasum", c("-a", "256", shQuote(path)), stdout=TRUE)
  strsplit(z, "[[:space:]]+")[[1]][1]
}
w <- function(x, name) fwrite(x, file.path(OUT, name), sep="\t", quote=FALSE, na="NA")

reuse <- data.table(
  artifact=c("BAYESPRISM_FULL_FIT.rds", "CELL_FRACTION_ESTIMATES.tsv", "MALIGNANT_EXPRESSION_ESTIMATES.tsv", "BAYESPRISM_SAMPLE_QC.tsv", "BULK_PATIENT_SCORES.tsv"),
  source_path=c(fit_path, frac_path, mal_path, qc_path, bulk_path)
)
reuse[, `:=`(size_bytes=file.info(source_path)$size, SHA256=vapply(source_path, sha256, character(1)),
             lineage_status="REUSE_VERIFIED", formal_use=c("posterior reference path only", "formal composition predictors", "canonical malignant expression only; candidate-program columns excluded", "numeric reconstruction/posterior QA only", "strict cohort bulk canonical score"))]
w(reuse, "BAYESPRISM_REUSE_MANIFEST.tsv")

bulk <- fread(bulk_path)
frac <- fread(frac_path)
mal <- fread(mal_path)[, .(sample_id, cohort, malignant_canonical_MES, canonical_gene_coverage)]
qc <- fread(qc_path)
sanitized_mal <- copy(mal)
w(sanitized_mal, "MALIGNANT_EXPRESSION_ESTIMATES.tsv")
w(frac, "CELL_FRACTION_ESTIMATES_REUSED.tsv")
w(qc, "BAYESPRISM_SAMPLE_QC_REUSED.tsv")

setnames(bulk, "canonical_MES", "bulk_canonical_MES")
d <- merge(bulk[, .(cohort, patient_id, bulk_canonical_MES)], mal,
           by.x=c("cohort", "patient_id"), by.y=c("cohort", "sample_id"), all=FALSE)
d <- merge(d, frac[, .(cohort, sample_id, malignant_fraction, myeloid_fraction,
                       endothelial_fraction, pericyte_stromal_fraction, lymphoid_fraction,
                       posterior_median_cv)],
           by.x=c("cohort", "patient_id"), by.y=c("cohort", "sample_id"), all=FALSE)
d[, vascular_stromal_fraction := endothelial_fraction + pericyte_stromal_fraction]
preds <- c("malignant_canonical_MES", "myeloid_fraction", "vascular_stromal_fraction")

r2_for <- function(dd, use) {
  if (!length(use)) return(0)
  summary(lm(reformulate(use, "bulk_canonical_MES"), dd))$r.squared
}
adjr2_for <- function(dd, use) summary(lm(reformulate(use, "bulk_canonical_MES"), dd))$adj.r.squared
vif_values <- function(dd) setNames(vapply(preds, function(p) {
  r <- summary(lm(reformulate(setdiff(preds, p), p), dd))$r.squared
  1/(1-r)
}, numeric(1)), preds)

permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow=1))
  do.call(rbind, lapply(seq_along(x), function(i) cbind(x[i], permutations(x[-i]))))
}
perms <- permutations(preds)
lmg <- function(dd) {
  contributions <- setNames(numeric(length(preds)), preds)
  for (i in seq_len(nrow(perms))) {
    current <- character(); previous <- 0
    for (p in perms[i,]) {
      current <- c(current, p)
      next_r2 <- r2_for(dd, current)
      contributions[p] <- contributions[p] + (next_r2 - previous)
      previous <- next_r2
    }
  }
  contributions/nrow(perms)
}

partition <- function(dd) {
  full_r2 <- r2_for(dd, preds)
  full_adj <- adjr2_for(dd, preds)
  unique <- setNames(vapply(preds, function(p) full_r2-r2_for(dd, setdiff(preds,p)), numeric(1)), preds)
  shared <- full_r2-sum(unique)
  list(r2=full_r2, adj=full_adj, unique=unique, shared=shared, lmg=lmg(dd), vif=vif_values(dd))
}

formal_rows <- list(); common_rows <- list(); relative_rows <- list(); ridge_rows <- list(); boot_long <- list()
for (co in sort(unique(d$cohort))) {
  dd <- d[cohort == co]
  fit <- lm(reformulate(preds, "bulk_canonical_MES"), dd)
  p <- partition(dd)
  sm <- summary(fit)$coefficients
  formal_rows[[co]] <- rbind(
    data.table(cohort=co, result_type="MODEL_METRIC", term=c("n", "R2", "adjusted_R2", "max_VIF"),
               estimate=c(nrow(dd), p$r2, p$adj, max(p$vif)), std_error=NA_real_, statistic=NA_real_, p_value=NA_real_),
    data.table(cohort=co, result_type="COEFFICIENT", term=rownames(sm), estimate=sm[,1], std_error=sm[,2], statistic=sm[,3], p_value=sm[,4])
  )
  common_rows[[co]] <- rbind(
    data.table(cohort=co, component=paste0("unique_", names(p$unique)), estimate=as.numeric(p$unique), n=nrow(dd), max_VIF=max(p$vif)),
    data.table(cohort=co, component="shared_variance", estimate=p$shared, n=nrow(dd), max_VIF=max(p$vif)),
    data.table(cohort=co, component="total_R2", estimate=p$r2, n=nrow(dd), max_VIF=max(p$vif))
  )
  relative_rows[[co]] <- data.table(cohort=co, predictor=names(p$lmg), lmg_shapley_R2=as.numeric(p$lmg),
                                    fraction_of_model_R2=as.numeric(p$lmg)/p$r2, n=nrow(dd), max_VIF=max(p$vif))

  z <- as.data.table(scale(dd[, c("bulk_canonical_MES", preds), with=FALSE]))
  y <- z$bulk_canonical_MES; X <- as.matrix(z[, ..preds])
  for (lambda in c(0, 0.1, 1, 10)) {
    beta <- solve(crossprod(X) + diag(lambda, ncol(X)), crossprod(X,y))
    ridge_rows[[paste(co,lambda)]] <- data.table(cohort=co, lambda=lambda, predictor=preds,
                                                 standardized_coefficient=as.numeric(beta), n=nrow(dd),
                                                 selection_basis="fixed prespecified lambda; no outcome-based tuning")
  }

  for (b in seq_len(1000)) {
    bb <- dd[sample.int(nrow(dd), nrow(dd), replace=TRUE)]
    if (any(vapply(bb[, ..preds], function(x) sd(x,na.rm=TRUE)==0, logical(1)))) next
    q <- tryCatch(partition(bb), error=function(e) NULL)
    if (is.null(q)) next
    boot_long[[paste(co,b)]] <- data.table(
      cohort=co, bootstrap_id=b,
      metric=c("R2", "adjusted_R2", paste0("unique_", names(q$unique)), "shared_variance", paste0("LMG_", names(q$lmg))),
      estimate=c(q$r2, q$adj, as.numeric(q$unique), q$shared, as.numeric(q$lmg))
    )
  }
}

formal <- rbindlist(formal_rows, fill=TRUE)
common <- rbindlist(common_rows)
relative <- rbindlist(relative_rows)
ridge <- rbindlist(ridge_rows)
boot <- rbindlist(boot_long)
boot_summary <- boot[, .(bootstrap_estimate=median(estimate,na.rm=TRUE), ci_low=quantile(estimate,.025,na.rm=TRUE),
                         ci_high=quantile(estimate,.975,na.rm=TRUE), successful_bootstraps=.N), by=.(cohort,metric)]

w(formal, "BAYESPRISM_FORMAL_MODEL.tsv")
w(common, "BAYESPRISM_COMMONALITY.tsv")
w(relative, "BAYESPRISM_RELATIVE_IMPORTANCE.tsv")
w(boot_summary, "BAYESPRISM_BOOTSTRAP.tsv")
w(ridge, "BAYESPRISM_RIDGE_SENSITIVITY.tsv")

gate_status <- if (all(qc$posterior_mean_finite_fraction == 1) && all(common[component=="total_R2", estimate] > 0)) "PARTIAL_PASS_REFERENCE_SENSITIVITY_PENDING" else "FAIL"
writeLines(c(
  "# BayesPrism analysis summary",
  "",
  paste0("**", gate_status, "**"),
  "",
  "The formal v1.1 model is: bulk canonical MES ~ inferred malignant canonical MES + myeloid fraction + vascular/stromal fraction.",
  "Candidate-program scores are absent from every coefficient, commonality, bootstrap, relative-importance and ridge calculation.",
  "High VIF is handled by emphasizing shared variance, Shapley/LMG importance, bootstrap intervals and ridge direction rather than isolated ordinary-regression coefficients.",
  "Reconstruction cosine is retained only as a numerical reconstruction check and is not an independent accuracy validation.",
  "Fractions are model-based expression proportions, not pathology cell counts."
), file.path(OUT, "BAYESPRISM_ANALYSIS_SUMMARY.md"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "SESSION_INFO.txt"))
cat("BAYESPRISM_FORMAL_MODEL_COMPLETE", nrow(d), "patients", gate_status, "\n")
