#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(readxl)
})

DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset = "data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset = "results")
REPO_ROOT <- Sys.getenv("GBM_MES_REPO_ROOT", unset = ".")
run_dir <- file.path(RESULT_ROOT, "clinical_context_reconciliation")
out_dir <- file.path(run_dir, "outputs")
log_dir <- file.path(run_dir, "logs")
work_dir <- file.path(run_dir, "work")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

score_file <- file.path(RESULT_ROOT, "clinical_icb", "02_full_icb_projection", "FULL_ICB_PATIENT_COMPONENT_SCORES.tsv")
cohort_file <- file.path(RESULT_ROOT, "clinical_icb", "01_full_icb_audit", "FULL_ICB_ANALYSIS_COHORT.tsv")
bench_file <- file.path(RESULT_ROOT, "clinical_specificity", "outputs", "P2_BENCHMARK_PATIENT_TABLE.tsv")
icb181_file <- file.path(DATA_ROOT, "clinical", "icb", "bulk_transcriptional_classifiers", "icb_cohort_unified.tsv")
counts_file <- file.path(DATA_ROOT, "clinical", "icb", "rnaseq_counts_tpms.xlsx")
sig_long_file <- file.path(REPO_ROOT, "resources", "signatures", "signature_gene_sets.tsv")
long_delta_file <- file.path(RESULT_ROOT, "clinical_icb", "03_longitudinal", "LONGITUDINAL_PATIENT_DELTAS.tsv")
long_gate_file <- file.path(RESULT_ROOT, "clinical_icb", "03_longitudinal", "LONGITUDINAL_LANDMARK_GATE.tsv")

write_tsv <- function(x, name) fwrite(as.data.table(x), file.path(out_dir, name), sep = "\t", quote = FALSE, na = "NA")
z <- function(x) as.numeric(scale(as.numeric(x)))
safe_num <- function(x) suppressWarnings(as.numeric(x))

cohort <- fread(cohort_file)
scores <- fread(score_file)
bench <- fread(bench_file)
dat <- merge(cohort, scores, by = intersect(c("participant_id", "sample_id"), intersect(names(cohort), names(scores))), all = FALSE)
dat <- merge(dat, bench, by = intersect(c("participant_id", "sample_id"), intersect(names(dat), names(bench))), all = FALSE, suffixes = c("", "_bench"))

pick <- function(candidates, nms = names(dat)) {
  hit <- candidates[candidates %in% nms]
  if (!length(hit)) stop("Missing required variable: ", paste(candidates, collapse = " / "))
  hit[[1]]
}

time_col <- pick(c("time_years", "OS_years", "os_years", "time"))
event_col <- pick(c("event", "OS_event", "os_event"))
mal_col <- pick(c("malignant_MES_sd", "malignant_mes_sd", "malignant_MES"))
eco_col <- pick(c("ecological_component_sd", "ecological_component"))
stable_col <- pick(c("stable_mixed_contribution_sd", "stable_mixed_contribution"))
bulk_col <- pick(c("bulk_MES_sd", "bulk_mes_sd", "bulk_MES"))

core <- data.table(
  participant_id = dat$participant_id,
  sample_id = dat$sample_id,
  time_years = safe_num(dat[[time_col]]),
  event = as.integer(dat[[event_col]]),
  malignant_MES = z(dat[[mal_col]]),
  ecological_MES = z(dat[[eco_col]]),
  stable_mixed = z(dat[[stable_col]]),
  bulk_MES = z(dat[[bulk_col]])
)

bench_map <- c(
  HLA_I = "HLA_I_ssGSEA",
  HLA_II = "HLA_II_ssGSEA",
  T_cell = "Neftel_T_cell_score",
  IFNg = "IFN_gamma_hallmark",
  Neftel_macrophage = "Neftel_macrophage_score",
  myeloid_fraction = "myeloid_fraction",
  hypoxia = "hypoxia_hallmark",
  matrix = "matrix_core_matrisome",
  vascular_stromal = "vascular_stromal_fraction",
  malignant_fraction = "malignant_fraction"
)
for (nm in names(bench_map)) {
  src <- bench_map[[nm]]
  core[[nm]] <- if (src %in% names(dat)) z(dat[[src]]) else NA_real_
}

## Outcome-independent new module implementation, frozen before outcome analysis.
raw <- as.data.table(read_excel(counts_file, sheet = "sTable6a_rnaseq_counts", skip = 1))
setnames(raw, 1:2, c("ensembl_id", "gene_symbol"))
count_cols <- setdiff(names(raw), c("ensembl_id", "gene_symbol"))
raw[, gene_symbol := toupper(trimws(as.character(gene_symbol)))]
raw <- raw[!is.na(gene_symbol) & gene_symbol != ""]
for (cc in count_cols) set(raw, j = cc, value = safe_num(raw[[cc]]))
gene_counts <- raw[, lapply(.SD, sum, na.rm = TRUE), by = gene_symbol, .SDcols = count_cols]
count_mat <- t(as.matrix(gene_counts[, ..count_cols]))
colnames(count_mat) <- gene_counts$gene_symbol
rownames(count_mat) <- count_cols
lib_size <- rowSums(count_mat, na.rm = TRUE)
log_cpm <- log2(sweep(count_mat, 1, lib_size / 1e6, "/") + 1)
gene_z <- scale(log_cpm)

sig_long <- fread(sig_long_file)
tam_genes <- unique(toupper(sig_long[signature_id == "MYELOID_ACTIVATED_TAM", gene_symbol]))
module_defs <- list(
  CYT_GZMA_PRF1 = c("GZMA", "PRF1"),
  cDC1_BATF3_XCR1_CLEC9A = c("BATF3", "XCR1", "CLEC9A"),
  MYELOID_ACTIVATED_TAM = tam_genes,
  OSM_OSMR_LIFR_STAT3 = c("OSM", "OSMR", "LIFR", "STAT3")
)

module_meta <- data.table(
  module = c("HLA_I", "HLA_II", "IFNg", "T_cell", "CYT_GZMA_PRF1", "cDC1_BATF3_XCR1_CLEC9A",
             "Neftel_macrophage", "myeloid_fraction", "MYELOID_ACTIVATED_TAM",
             "hypoxia", "matrix", "vascular_stromal", "malignant_fraction", "OSM_OSMR_LIFR_STAT3"),
  family = c(rep("IMMUNE_COMPETENT", 6), rep("SUPPRESSIVE_STRUCTURAL", 7), "CONTEXTUAL_MECHANISTIC"),
  priority = c(rep("PRIMARY", 4), rep("SECONDARY", 2), rep("SECONDARY", 7), "CONTEXTUAL"),
  source = c("Frozen project HLA-I ssGSEA", "Frozen project HLA-II ssGSEA",
             "Frozen project Hallmark IFN-gamma", "Frozen project Neftel T-cell",
             "Rooney et al. Cell 2015", "Published cDC1 marker set",
             "Frozen project Neftel macrophage", "Frozen BayesPrism myeloid fraction",
             "Frozen project MYELOID_ACTIVATED_TAM", "Frozen project Hallmark hypoxia",
             "Frozen project core matrisome", "Frozen BayesPrism vascular/stromal fraction",
             "Frozen BayesPrism malignant fraction", "Hara et al. Cancer Cell 2021 contextual axis")
)
module_meta[, outcome_independent_inclusion := "YES"]
module_meta[, score_definition := fifelse(module %in% names(module_defs),
  "Equal-weight mean of gene-wise z scores on log2(CPM+1), then patient z",
  "Existing frozen project score; patient z-standardized")]
module_meta[, genes_original := NA_integer_]
module_meta[, genes_mapped := NA_integer_]
module_meta[, mapping_rate := NA_real_]
module_meta[, genes_requested := NA_character_]
module_meta[, genes_available := NA_character_]

for (mn in names(module_defs)) {
  req <- unique(module_defs[[mn]])
  hit <- intersect(req, colnames(gene_z))
  module_meta[module == mn, `:=`(
    genes_original = length(req), genes_mapped = length(hit),
    mapping_rate = length(hit) / length(req),
    genes_requested = paste(req, collapse = ";"),
    genes_available = paste(hit, collapse = ";")
  )]
  if (length(hit) / length(req) >= 0.80) {
    sc <- rowMeans(gene_z[, hit, drop = FALSE], na.rm = TRUE)
    mapped <- data.table(sample_id = rownames(gene_z), value = z(sc))
    setnames(mapped, "value", mn)
    core <- merge(core, mapped, by = "sample_id", all.x = TRUE)
  } else {
    core[[mn]] <- NA_real_
  }
}
module_meta[, certification_status := fifelse(
  is.na(mapping_rate) | mapping_rate >= 0.80, "CERTIFIED_FROZEN_OR_PUBLISHED", "NOT_EVALUABLE_MAPPING_BELOW_80_PERCENT")]
write_tsv(module_meta, "P2_FUNCTIONAL_ECOLOGY_MODULE_REGISTRY.tsv")

stopifnot(nrow(core) == 40L, sum(core$event) == 38L)

## P1 variable availability and Ghannam181 non-overlap feasibility audit.
icb181 <- fread(icb181_file)
count_nonmissing <- function(cols) {
  hit <- intersect(cols, names(icb181))
  if (!length(hit)) return(0L)
  max(vapply(hit, function(cc) sum(!is.na(icb181[[cc]]) & as.character(icb181[[cc]]) != ""), integer(1)))
}
availability <- data.table(
  variable = c("participant_id", "OS_time", "OS_event", "pretreatment_RNA_sample_id",
               "bulk_MES_subtype", "bulk_MES_continuous", "direct_Ghannam_malignant_MES_like",
               "malignant_MES_source_score", "ecological_MES_source_score", "stable_mixed",
               "HLA_I", "HLA_II", "T_cell", "IFNg", "immune_infiltration", "source_study_membership",
               "sample_timing", "nonoverlap_complete_pretreatment_RNA_outcome"),
  availability = c("DIRECT", "DIRECT", "DIRECT", "DIRECT",
                   "DIRECT", "RECONSTRUCTABLE_FROZEN", "PUBLICATION_ONLY",
                   "RECONSTRUCTABLE_FROZEN", "RECONSTRUCTABLE_FROZEN", "RECONSTRUCTABLE_FROZEN",
                   "DIRECT", "DIRECT", "RECONSTRUCTABLE_FROZEN", "RECONSTRUCTABLE_FROZEN",
                   "RECONSTRUCTABLE_FROZEN", "NOT_AVAILABLE", "DIRECT", "NOT_AVAILABLE"),
  evaluable_n = c(nrow(icb181),
                  count_nonmissing(c("osicb", "os", "OS", "overall_survival", "os_months", "time_years")),
                  count_nonmissing(c("Deceased", "event", "OS_event", "status")),
                  count_nonmissing(c("sample_id_legacy_pre", "pre_rna_id", "pre_sample_id", "sample_id_pre")),
                  count_nonmissing(c("verhaak_subtype_pre", "VERHAAK_GLIOBLASTOMA_MESENCHYMAL_pre", "pre_verhaak", "verhaak", "MES_status")),
                  40L, 0L, 40L, 40L, 40L,
                  count_nonmissing(c("GOCC_MHC_CLASS_I_PROTEIN_COMPLEX_pre", "pre_hla_i", "HLA_I")), count_nonmissing(c("GOCC_MHC_CLASS_II_PROTEIN_COMPLEX_pre", "pre_hla_ii", "HLA_II")),
                  40L, 40L, 40L, 0L,
                  count_nonmissing(c("collection_date_icb_pre", "collection_date_icb_post", "sample_id_legacy_pre", "sample_id_legacy_post", "sample_timing")), 0L),
  role = c("identity", "endpoint", "endpoint", "sample linkage", "published classifier",
           "frozen reconstruction", "reconciliation target", "frozen source decomposition",
           "frozen source decomposition", "frozen source decomposition", "functional benchmark",
           "functional benchmark", "functional benchmark", "functional benchmark",
           "functional benchmark", "overlap audit", "time-zero audit", "external transport"),
  note = c(
    "The aggregated Ghannam reference contains 181 participants.",
    "Patient-level survival is available in the aggregated reference.",
    "Patient-level event indicator is available in the aggregated reference.",
    "Pretreatment RNA identifiers are available only for the formal 40-patient continuous-score cohort.",
    "Pretreatment classifier annotation is available only where pretreatment RNA exists.",
    "Frozen project score exists for formal ICB40.",
    "No distinct patient-level Ghannam malignant-MES-like variable certified in the 181 annotation.",
    "Frozen source-aware score exists for formal ICB40.",
    "Frozen source-aware score exists for formal ICB40.",
    "Frozen source-aware score exists for formal ICB40.",
    "Direct project annotation on formal ICB40.",
    "Direct project annotation on formal ICB40.",
    "Frozen outcome-independent project module on ICB40.",
    "Frozen Hallmark module on ICB40.",
    "Frozen deconvolution/module measures on ICB40.",
    "No certified source-study membership field enabling independence adjudication.",
    "Pre/post identifiers exist only for a subset and do not certify a nonoverlap external cohort.",
    "181 minus 40 is not assumed to be an analyzable independent cohort."
  )
)
write_tsv(availability, "P1_GHANNAM_VARIABLE_AVAILABILITY.tsv")

## Cox helpers.
ph_extract <- function(fit) {
  zz <- tryCatch(cox.zph(fit, transform = "km"), error = function(e) NULL)
  if (is.null(zz)) return(list(global = NA_real_, terms = numeric()))
  tab <- as.data.frame(zz$table)
  list(global = if ("GLOBAL" %in% rownames(tab)) tab["GLOBAL", "p"] else NA_real_,
       terms = setNames(tab[setdiff(rownames(tab), "GLOBAL"), "p"], setdiff(rownames(tab), "GLOBAL")))
}

max_dfbeta <- function(fit) {
  rr <- tryCatch(residuals(fit, type = "dfbeta"), error = function(e) NULL)
  if (is.null(rr)) return(NA_real_)
  max(abs(rr), na.rm = TRUE)
}

design_diagnostics <- function(d, vars) {
  x <- as.matrix(d[, ..vars])
  x <- x[complete.cases(x), , drop = FALSE]
  if (!nrow(x) || ncol(x) < 2) return(list(max_vif = NA_real_, condition_number = NA_real_))
  cm <- cor(x)
  vif <- tryCatch(diag(solve(cm)), error = function(e) rep(Inf, ncol(cm)))
  list(max_vif = max(vif, na.rm = TRUE), condition_number = kappa(scale(x), exact = TRUE))
}

fit_model <- function(d, model_id, vars, benchmark = NA_character_) {
  cc <- d[complete.cases(d[, c("time_years", "event", vars), with = FALSE])]
  f <- as.formula(paste("Surv(time_years, event) ~", paste(vars, collapse = " + ")))
  fit <- coxph(f, data = cc, x = TRUE, model = TRUE, singular.ok = TRUE)
  sm <- summary(fit)
  ph <- ph_extract(fit)
  dd <- design_diagnostics(cc, vars)
  rows <- rbindlist(lapply(seq_along(coef(fit)), function(j) {
    tr <- names(coef(fit))[j]
    data.table(
      model_id = model_id, benchmark = benchmark, n = nrow(cc), events = sum(cc$event),
      term = tr, beta = unname(coef(fit)[j]), HR = exp(unname(coef(fit)[j])),
      CI_low = sm$conf.int[j, "lower .95"], CI_high = sm$conf.int[j, "upper .95"],
      Wald_P = sm$coefficients[j, "Pr(>|z|)"],
      PH_global_P = ph$global, PH_term_P = unname(ph$terms[tr]),
      max_abs_DFBETA = max_dfbeta(fit), max_VIF = dd$max_vif,
      condition_number = dd$condition_number,
      multiplicity_family = "MECHANISTIC_INFORMATION_ATTRIBUTION",
      standardized_unit = "PER_1_WITHIN_COHORT_SD"
    )
  }))
  list(fit = fit, data = cc, rows = rows)
}

p1_specs <- list(
  M0 = c("malignant_MES"),
  M1 = c("ecological_MES"),
  M2 = c("malignant_MES", "ecological_MES"),
  M3a_HLA_I = c("malignant_MES", "ecological_MES", "HLA_I"),
  M3b_HLA_II = c("malignant_MES", "ecological_MES", "HLA_II"),
  M3c_T_cell = c("malignant_MES", "ecological_MES", "T_cell"),
  M3d_IFNg = c("malignant_MES", "ecological_MES", "IFNg")
)
p1_models <- lapply(names(p1_specs), function(mid) fit_model(core, mid, p1_specs[[mid]],
  if (length(p1_specs[[mid]]) == 3) p1_specs[[mid]][3] else NA_character_))
names(p1_models) <- names(p1_specs)
p1_model_rows <- rbindlist(lapply(p1_models, `[[`, "rows"), fill = TRUE)

tv_estimates <- function(d, model_id, benchmark = NA_character_) {
  vars <- c("malignant_MES", "ecological_MES", if (!is.na(benchmark)) benchmark)
  cc <- d[complete.cases(d[, c("time_years", "event", vars), with = FALSE])]
  rhs <- paste(c("malignant_MES", if (!is.na(benchmark)) benchmark, "tt(ecological_MES)"), collapse = " + ")
  fit <- tryCatch(coxph(as.formula(paste("Surv(time_years,event) ~", rhs)), data = cc,
                        tt = function(x, t, ...) cbind(x, x * log(pmax(t, 1e-6) / (6/12))),
                        x = TRUE), error = function(e) NULL)
  if (is.null(fit)) return(data.table())
  cf <- coef(fit)
  vc <- vcov(fit)
  eco_idx <- grep("tt\\(ecological_MES\\)", names(cf))
  if (length(eco_idx) < 2) return(data.table())
  rbindlist(lapply(c(6, 12, 18), function(mo) {
    l <- log(mo / 6)
    w <- rep(0, length(cf)); w[eco_idx[1]] <- 1; w[eco_idx[2]] <- l
    b <- sum(w * cf); se <- sqrt(drop(t(w) %*% vc %*% w))
    data.table(model_id = model_id, benchmark = benchmark, month = mo,
               ecological_beta_t = b, ecological_HR_t = exp(b),
               CI_low = exp(b - 1.96 * se), CI_high = exp(b + 1.96 * se),
               ecological_favorable_at_time = b < 0)
  }))
}

tv_rows <- rbindlist(lapply(c("M2", grep("^M3", names(p1_specs), value = TRUE)), function(mid) {
  bm <- if (length(p1_specs[[mid]]) == 3) p1_specs[[mid]][3] else NA_character_
  tv_estimates(core, mid, bm)
}), fill = TRUE)

set.seed(2026081201)
B <- 1000L
p1_boot <- vector("list", B * 5L)
kk <- 1L
boot_model_ids <- c("M2", grep("^M3", names(p1_specs), value = TRUE))
for (b in seq_len(B)) {
  ids <- sample(seq_len(nrow(core)), replace = TRUE)
  db <- core[ids]
  for (mid in boot_model_ids) {
    vars <- p1_specs[[mid]]
    cc <- db[complete.cases(db[, c("time_years", "event", vars), with = FALSE])]
    fit <- tryCatch(coxph(as.formula(paste("Surv(time_years,event) ~", paste(vars, collapse = " + "))),
                          data = cc, ties = "efron", singular.ok = TRUE), error = function(e) NULL)
    ok <- !is.null(fit) && all(c("malignant_MES", "ecological_MES") %in% names(coef(fit))) &&
      all(is.finite(coef(fit)[c("malignant_MES", "ecological_MES")]))
    if (ok) {
      cf <- coef(fit)
      ph <- ph_extract(fit)
      bm <- if (length(vars) == 3) vars[3] else NA_character_
      p1_boot[[kk]] <- data.table(
        iteration = b, model_id = mid, benchmark = bm, estimable = TRUE,
        beta_malignant = unname(cf["malignant_MES"]),
        beta_ecological = unname(cf["ecological_MES"]),
        beta_benchmark = if (!is.na(bm) && bm %in% names(cf)) unname(cf[bm]) else NA_real_,
        malignant_favorable = cf["malignant_MES"] < 0,
        ecological_favorable = cf["ecological_MES"] < 0,
        ecological_more_favorable_than_malignant = cf["ecological_MES"] < cf["malignant_MES"],
        PH_compatible = is.na(ph$global) | ph$global >= 0.05,
        seed = 2026081201L
      )
    } else {
      p1_boot[[kk]] <- data.table(iteration = b, model_id = mid,
        benchmark = if (length(vars) == 3) vars[3] else NA_character_, estimable = FALSE,
        beta_malignant = NA_real_, beta_ecological = NA_real_, beta_benchmark = NA_real_,
        malignant_favorable = NA, ecological_favorable = NA,
        ecological_more_favorable_than_malignant = NA, PH_compatible = NA,
        seed = 2026081201L)
    }
    kk <- kk + 1L
  }
}
p1_boot <- rbindlist(p1_boot, fill = TRUE)
write_tsv(p1_boot, "P1_CONDITIONAL_ATTRIBUTION_BOOTSTRAP.tsv")

m2_beta_m <- p1_model_rows[model_id == "M2" & term == "malignant_MES", beta][1]
m2_beta_e <- p1_model_rows[model_id == "M2" & term == "ecological_MES", beta][1]
retention <- rbindlist(lapply(boot_model_ids, function(mid) {
  mr <- p1_model_rows[model_id == mid & term == "malignant_MES"]
  er <- p1_model_rows[model_id == mid & term == "ecological_MES"]
  bs <- p1_boot[model_id == mid & estimable == TRUE]
  tv <- tv_rows[model_id == mid]
  data.table(
    model_id = mid, benchmark = unique(er$benchmark)[1], n = er$n[1], events = er$events[1],
    beta_malignant = mr$beta[1], beta_ecological = er$beta[1],
    malignant_attenuation_vs_M2 = 1 - mr$beta[1] / m2_beta_m,
    ecological_attenuation_vs_M2 = 1 - er$beta[1] / m2_beta_e,
    malignant_sign_retained = sign(mr$beta[1]) == sign(m2_beta_m),
    ecological_sign_retained = sign(er$beta[1]) == sign(m2_beta_e),
    bootstrap_estimable_fraction = mean(p1_boot[model_id == mid]$estimable),
    bootstrap_malignant_favorable_fraction = mean(bs$malignant_favorable, na.rm = TRUE),
    bootstrap_ecological_favorable_fraction = mean(bs$ecological_favorable, na.rm = TRUE),
    bootstrap_ecological_more_favorable_fraction = mean(bs$ecological_more_favorable_than_malignant, na.rm = TRUE),
    bootstrap_PH_compatible_fraction = mean(bs$PH_compatible, na.rm = TRUE),
    ecological_beta_6m = tv[month == 6, ecological_beta_t][1],
    ecological_beta_12m = tv[month == 12, ecological_beta_t][1],
    ecological_beta_18m = tv[month == 18, ecological_beta_t][1],
    no_fixed_time_direction_reversal_6_18m = all(tv$ecological_beta_t < 0, na.rm = TRUE)
  )
}))
write_tsv(p1_model_rows, "P1_CONDITIONAL_ATTRIBUTION_MODELS.tsv")
write_tsv(retention, "P1_SOURCE_INFORMATION_RETENTION.tsv")

## P2 source-aware functional attenuation and paired bootstrap.
eval_modules <- module_meta[certification_status == "CERTIFIED_FROZEN_OR_PUBLISHED", module]
p2_models <- rbindlist(lapply(eval_modules, function(bm) {
  fit_model(core, paste0("M2_PLUS_", bm), c("malignant_MES", "ecological_MES", bm), bm)$rows
}), fill = TRUE)

set.seed(2026081202)
p2_boot_list <- vector("list", B * length(eval_modules))
kk <- 1L
for (b in seq_len(B)) {
  ids <- sample(seq_len(nrow(core)), replace = TRUE)
  db <- core[ids]
  for (bm in eval_modules) {
    vars <- c("malignant_MES", "ecological_MES", bm)
    cc <- db[complete.cases(db[, c("time_years", "event", vars), with = FALSE])]
    fit <- tryCatch(coxph(as.formula(paste("Surv(time_years,event) ~", paste(vars, collapse = " + "))),
                          data = cc, ties = "efron", singular.ok = TRUE), error = function(e) NULL)
    ok <- !is.null(fit) && all(vars %in% names(coef(fit))) && all(is.finite(coef(fit)[vars]))
    if (ok) {
      cf <- coef(fit); ph <- ph_extract(fit)
      p2_boot_list[[kk]] <- data.table(
        iteration = b, benchmark = bm, estimable = TRUE,
        beta_malignant = unname(cf["malignant_MES"]),
        beta_ecological = unname(cf["ecological_MES"]),
        beta_benchmark = unname(cf[bm]),
        malignant_favorable = cf["malignant_MES"] < 0,
        ecological_favorable = cf["ecological_MES"] < 0,
        ecological_more_favorable_than_malignant = cf["ecological_MES"] < cf["malignant_MES"],
        PH_compatible = is.na(ph$global) | ph$global >= 0.05,
        seed = 2026081202L)
    } else {
      p2_boot_list[[kk]] <- data.table(iteration = b, benchmark = bm, estimable = FALSE,
        beta_malignant = NA_real_, beta_ecological = NA_real_, beta_benchmark = NA_real_,
        malignant_favorable = NA, ecological_favorable = NA,
        ecological_more_favorable_than_malignant = NA, PH_compatible = NA, seed = 2026081202L)
    }
    kk <- kk + 1L
  }
}
p2_boot <- rbindlist(p2_boot_list, fill = TRUE)
write_tsv(p2_boot, "P2_SOURCE_AWARE_FUNCTIONAL_BOOTSTRAP.tsv")

p2_att <- rbindlist(lapply(eval_modules, function(bm) {
  rows <- p2_models[benchmark == bm]
  mr <- rows[term == "malignant_MES"]; er <- rows[term == "ecological_MES"]; br <- rows[term == bm]
  bs <- p2_boot[benchmark == bm & estimable == TRUE]
  meta <- module_meta[module == bm]
  data.table(
    benchmark = bm, family = meta$family[1], priority = meta$priority[1],
    n = er$n[1], events = er$events[1],
    beta_malignant = mr$beta[1], HR_malignant = mr$HR[1], malignant_P = mr$Wald_P[1],
    beta_ecological = er$beta[1], HR_ecological = er$HR[1], ecological_P = er$Wald_P[1],
    beta_benchmark = br$beta[1], HR_benchmark = br$HR[1], benchmark_P = br$Wald_P[1],
    malignant_attenuation_vs_M2 = 1 - mr$beta[1] / m2_beta_m,
    ecological_attenuation_vs_M2 = 1 - er$beta[1] / m2_beta_e,
    ecological_sign_retained = sign(er$beta[1]) == sign(m2_beta_e),
    bootstrap_estimable_fraction = mean(p2_boot[benchmark == bm]$estimable),
    bootstrap_ecological_favorable_fraction = mean(bs$ecological_favorable, na.rm = TRUE),
    bootstrap_malignant_favorable_fraction = mean(bs$malignant_favorable, na.rm = TRUE),
    bootstrap_ecological_more_favorable_fraction = mean(bs$ecological_more_favorable_than_malignant, na.rm = TRUE),
    PH_global_P = er$PH_global_P[1], PH_ecological_term_P = er$PH_term_P[1],
    max_abs_DFBETA = er$max_abs_DFBETA[1], max_VIF = er$max_VIF[1],
    condition_number = er$condition_number[1],
    nonseparable_due_to_collinearity = er$max_VIF[1] >= 5 | er$condition_number[1] >= 30
  )
}))
p2_att[, benchmark_BH_FDR := p.adjust(benchmark_P, method = "BH")]
write_tsv(p2_att, "P2_SOURCE_AWARE_FUNCTIONAL_ATTENUATION.tsv")

## Outcome-independent correlation architecture.
corr_vars <- unique(c("malignant_MES", "ecological_MES", "stable_mixed", eval_modules))
corr_vars <- corr_vars[corr_vars %in% names(core)]
set.seed(2026081203)
corr_rows <- list(); kk <- 1L
for (i in seq_along(corr_vars)) for (j in i:length(corr_vars)) {
  a <- corr_vars[i]; b <- corr_vars[j]
  cc <- core[complete.cases(core[, c(a, b), with = FALSE])]
  rho <- suppressWarnings(cor(cc[[a]], cc[[b]], method = "spearman"))
  pv <- if (a == b) 0 else suppressWarnings(cor.test(cc[[a]], cc[[b]], method = "spearman", exact = FALSE)$p.value)
  boots <- replicate(B, {
    id <- sample(seq_len(nrow(cc)), replace = TRUE)
    suppressWarnings(cor(cc[[a]][id], cc[[b]][id], method = "spearman"))
  })
  corr_rows[[kk]] <- data.table(variable_1 = a, variable_2 = b, n = nrow(cc), rho = rho, P = pv,
    bootstrap_CI_low = quantile(boots, .025, na.rm = TRUE),
    bootstrap_CI_high = quantile(boots, .975, na.rm = TRUE),
    bootstrap_positive_fraction = mean(boots > 0, na.rm = TRUE), seed = 2026081203L)
  kk <- kk + 1L
}
corr_dt <- rbindlist(corr_rows)
corr_dt[, BH_FDR := p.adjust(P, method = "BH")]
write_tsv(corr_dt, "P2_FUNCTIONAL_ECOLOGY_CORRELATION_MATRIX.tsv")

partial_spearman <- function(x, y, control) {
  rx <- residuals(lm(rank(x, na.last = "keep") ~ rank(control, na.last = "keep")))
  ry <- residuals(lm(rank(y, na.last = "keep") ~ rank(control, na.last = "keep")))
  suppressWarnings(cor(rx, ry))
}
partial_rows <- list(); kk <- 1L
for (bm in eval_modules) {
  for (target in c("ecological_MES", "malignant_MES")) {
    control <- if (target == "ecological_MES") "malignant_MES" else "ecological_MES"
    cc <- core[complete.cases(core[, c(target, bm, control), with = FALSE])]
    rho <- partial_spearman(cc[[target]], cc[[bm]], cc[[control]])
    boots <- replicate(B, {
      id <- sample(seq_len(nrow(cc)), replace = TRUE)
      tryCatch(partial_spearman(cc[[target]][id], cc[[bm]][id], cc[[control]][id]), error = function(e) NA_real_)
    })
    partial_rows[[kk]] <- data.table(
      target = target, benchmark = bm, conditioned_on = control, n = nrow(cc),
      partial_spearman_rho = rho,
      bootstrap_CI_low = quantile(boots, .025, na.rm = TRUE),
      bootstrap_CI_high = quantile(boots, .975, na.rm = TRUE),
      bootstrap_positive_fraction = mean(boots > 0, na.rm = TRUE),
      seed = 2026081203L)
    kk <- kk + 1L
  }
}
partial_dt <- rbindlist(partial_rows)
write_tsv(partial_dt, "P2_PARTIAL_CORRELATION_ARCHITECTURE.tsv")

## P3 continuous longitudinal source geometry.
ld <- fread(long_delta_file)
ld_pick <- function(candidates) {
  hit <- candidates[candidates %in% names(ld)]
  if (!length(hit)) stop("Longitudinal variable missing: ", paste(candidates, collapse = " / "))
  hit[[1]]
}
idc <- ld_pick(c("participant_id", "patient_id"))
mpre <- ld_pick(c("malignant_MES_z_pre", "malignant_pre", "pre_malignant_MES", "malignant_MES_pre"))
mpost <- ld_pick(c("malignant_MES_z_post", "malignant_post", "post_malignant_MES", "malignant_MES_post"))
epre <- ld_pick(c("ecological_component_z_pre", "ecological_pre", "pre_ecological_component", "ecological_component_pre"))
epost <- ld_pick(c("ecological_component_z_post", "ecological_post", "post_ecological_component", "ecological_component_post"))
traj <- data.table(
  participant_id = ld[[idc]],
  malignant_pre = safe_num(ld[[mpre]]), ecological_pre = safe_num(ld[[epre]]),
  malignant_post = safe_num(ld[[mpost]]), ecological_post = safe_num(ld[[epost]])
)
traj[, `:=`(
  delta_malignant = malignant_post - malignant_pre,
  delta_ecological = ecological_post - ecological_pre,
  baseline_coupling_position = (malignant_pre + ecological_pre) / sqrt(2),
  baseline_rewiring_position = (malignant_pre - ecological_pre) / sqrt(2),
  post_coupling_position = (malignant_post + ecological_post) / sqrt(2),
  post_rewiring_position = (malignant_post - ecological_post) / sqrt(2)
)]
traj[, `:=`(
  delta_coupled = (delta_malignant + delta_ecological) / sqrt(2),
  delta_rewiring = (delta_malignant - delta_ecological) / sqrt(2),
  vector_magnitude = sqrt(delta_malignant^2 + delta_ecological^2),
  vector_angle_radians = atan2(delta_ecological, delta_malignant),
  vector_angle_degrees = atan2(delta_ecological, delta_malignant) * 180 / pi
)]
write_tsv(traj, "P3_LONGITUDINAL_SOURCE_TRAJECTORIES.tsv")

geom <- rbindlist(list(
  cbind(data.table(record_type = "PATIENT"), traj),
  data.table(record_type = "SUMMARY", participant_id = "ALL_19",
    malignant_pre = median(traj$malignant_pre), ecological_pre = median(traj$ecological_pre),
    malignant_post = median(traj$malignant_post), ecological_post = median(traj$ecological_post),
    delta_malignant = median(traj$delta_malignant), delta_ecological = median(traj$delta_ecological),
    baseline_coupling_position = median(traj$baseline_coupling_position),
    baseline_rewiring_position = median(traj$baseline_rewiring_position),
    post_coupling_position = median(traj$post_coupling_position),
    post_rewiring_position = median(traj$post_rewiring_position),
    delta_coupled = median(traj$delta_coupled), delta_rewiring = median(traj$delta_rewiring),
    vector_magnitude = median(traj$vector_magnitude),
    vector_angle_radians = NA_real_, vector_angle_degrees = NA_real_)
), fill = TRUE)
write_tsv(geom, "P3_COUPLING_REWIRING_GEOMETRY.tsv")

gate <- fread(long_gate_file)
outcome_gate <- data.table(
  analysis = c("delta_coupled_vs_post_sample_outcome", "delta_rewiring_vs_post_sample_outcome"),
  status = "NOT_EVALUABLE",
  reason = "MECHANISTIC_EXPLORATORY_ONLY; frozen landmark gate is NOT_MODELED_UNDERPOWERED_AND_SELECTED_BY_REOPERATION",
  nominal_pairs = nrow(traj), landmark_evaluable_n = 18L, landmark_events = 17L,
  model_run = "NO", cutoff_used = "NO", reference_source = long_gate_file
)
write_tsv(outcome_gate, "P3_LONGITUDINAL_OUTCOME_EXPLORATORY.tsv")

## Machine-readable run facts used by integration.
facts <- list(
  n_icb40 = nrow(core), events_icb40 = sum(core$event),
  m2_beta_malignant = m2_beta_m, m2_beta_ecological = m2_beta_e,
  m2_boot_eco_fav = retention[model_id == "M2", bootstrap_ecological_favorable_fraction],
  m2_boot_eco_more = retention[model_id == "M2", bootstrap_ecological_more_favorable_fraction],
  m3_models_meeting_gate = sum(retention[grepl("^M3", model_id),
    bootstrap_ecological_favorable_fraction >= 0.80 &
    bootstrap_ecological_more_favorable_fraction >= 0.80 &
    no_fixed_time_direction_reversal_6_18m], na.rm = TRUE),
  evaluable_modules = length(eval_modules),
  longitudinal_pairs = nrow(traj)
)
saveRDS(facts, file.path(work_dir, "r7_2_primary_facts.rds"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "R7_2_PRIMARY_SESSION_INFO.txt"))
writeLines("R7.2 primary analysis completed successfully.", file.path(log_dir, "R7_2_PRIMARY_COMPLETION.log"))
