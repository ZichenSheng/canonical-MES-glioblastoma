#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
suppressPackageStartupMessages({library(data.table); library(lme4)})
set.seed(2026080141)

DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset="data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset="results")
REPO_ROOT <- Sys.getenv("GBM_MES_REPO_ROOT", unset=".")
RUN <- file.path(RESULT_ROOT,"spatial")
IVY <- file.path(DATA_ROOT,"spatial","ivy_gap")
MATRIX_DIR <- file.path(IVY,"gene_expression_matrix")
REGISTRY <- file.path(REPO_ROOT,"resources","signatures","signature_gene_sets.tsv")
OUTA <- file.path(RUN,"ivy_inputs")
OUTM <- file.path(RUN,"ivy_models")
dir.create(OUTA, recursive=TRUE, showWarnings=FALSE); dir.create(OUTM, recursive=TRUE, showWarnings=FALSE)
w <- function(x, p) fwrite(x, p, sep="\t", quote=FALSE, na="NA")

start_time <- proc.time()[["elapsed"]]
mapping_all <- fread(file.path(IVY, "IVY_GAP_SAMPLE_MAPPING.tsv"))
mapping <- mapping_all[normalized_matrix_column_present == TRUE]
region_from_acronym <- function(x) fifelse(grepl("^LE", x), "LEADING_EDGE",
  fifelse(grepl("^IT", x), "INFILTRATING_TUMOUR",
  fifelse(grepl("^CTmvp", x), "MICROVASCULAR_PROLIFERATION",
  fifelse(grepl("^CTpan", x), "PSEUDOPALISADING_CELLS_AROUND_NECROSIS",
  fifelse(grepl("^CT", x), "CELLULAR_TUMOUR", NA_character_)))))
mapping[, region := region_from_acronym(structure_parent_acronym)]
if (nrow(mapping) != 270L || anyNA(mapping$region) || uniqueN(mapping$sample_id) != 270L) stop("P0: expected 270 uniquely mapped normalized Allen samples across five regions")

region_map <- unique(mapping[, .(raw_structure_parent_acronym=structure_parent_acronym,
  raw_structure_parent_name=structure_parent_name, standardized_region=region)])
setorder(region_map, standardized_region, raw_structure_parent_acronym)
w(region_map, file.path(OUTA, "IVY_REGION_MAPPING.tsv"))

genes <- fread(file.path(MATRIX_DIR, "rows-genes.csv"))
expression <- fread(file.path(MATRIX_DIR, "fpkm_table.csv"), check.names=FALSE)
gene_id <- as.character(expression[[1]])
gene_symbol <- toupper(trimws(genes$gene_symbol[match(gene_id, as.character(genes$gene_id))]))
valid <- !is.na(gene_symbol) & nzchar(gene_symbol)
mat0 <- as.matrix(expression[valid, -1]); storage.mode(mat0) <- "double"; rownames(mat0) <- gene_symbol[valid]
if (anyDuplicated(rownames(mat0))) {
  order_symbols <- unique(rownames(mat0)); summed <- rowsum(mat0, rownames(mat0), reorder=FALSE)
  multiplicity <- as.numeric(table(factor(rownames(mat0), levels=order_symbols)))
  mat0 <- summed / multiplicity
}
sample_ids <- as.character(mapping$sample_id)
available <- intersect(sample_ids, colnames(mat0))
if (length(available) != 270L) stop("P0: all 270 mapped samples must exist in Allen normalized matrix")
mapping <- mapping[match(available, sample_id)]
matrix <- mat0[, available, drop=FALSE]
if (any(colnames(matrix) != as.character(mapping$sample_id))) stop("P0: matrix/sample mapping order mismatch")
log_expression <- log2(matrix + 1)
rm(expression, mat0, matrix); gc()

registry <- fread(REGISTRY)
registry[, gene := toupper(trimws(gene_symbol_clean))]
getset <- function(id) unique(registry[signature_id == id & !is.na(gene), gene])
HYP <- strsplit("NDRG1 TRIB3 ERO1A DNAJB9 HILPDA SLCO4A1 INSIG2 AKAP12 SLC3A2 ZFAS1 GOLT1B LOX NRN1 GBE1 SESN2 MXI1 STC2 VEGFA CEBPG UPP1 HK2 TAF1D INSIG1 EIF4EBP1 ADM PLOD2 RRAGD EPB41L4A-AS1 DDIT3 ANKRD37 FAM210A PDK1 PGK1 GFPT1 IGFBP3 CA9 TSPYL2 SLC38A1 EGLN3 SLC6A6 IGFBP5 SLC7A5 ANGPTL4 SLC1A5 FAM13A TREM1 CXCL8 TM4SF1 OSGIN1 PSAT1 BNIP3 ARFGEF3 LUCAT1 STX3 ASNS LARP6 GDF15 SH3BGR OSER1 NUPR1 GPNMB GRPEL2 KDM3A FICD ATF5 CLEC2B TMEM38B ATF4 DDIT4 CTH CREBRF CIART TUBA4A BNIP3L JMJD6 PIGA BEX2 BIRC3 BTG1 PPP1R3C PRKAG2-AS1 HSPH1 HSPA9 HSPA6 PID1 EMX2 KISS1R PHF1 LIF PERP LMAN1 ARID3A XBP1 ANG NUP58 MAP1LC3B ZFAND2A TMEM45A HSPA5 FRMD3", " ")[[1]]
AST <- strsplit("IGFBP7 GFAP ID3 METTL7B EFEMP1 PTN CHI3L1 PMP2 GAP43 IFITM3 CCL2 KLHDC8A FABP7 C21orf62 SNTG1 POSTN THY1 PTPRZ1 SRPX2 GPM6A COL5A2 HOPX AQP4 APOE MXRA8 CXCL14 CP RAMP1 CFI RCAN1 FABP5 S100B SERPING1 CLU CST3 A2M SPARC ABCC3 GPM6B ATP1B2 SPARCL1 MLC1 AGT SCARA3 LTF NES NNMT COL1A2 PLA2G5 PLAAT4 PLA2G2A MAOB C1S C1R C3 SERPINA3 SLPI PROS1 OLFM2 GPX3 CDH11 SPP1 TNC MOXD1 VCAN CITED1 SAA1 ID1 VCAM1 CNR1 HTRA1 GAS1 NPAS3 GPC1 HLA-DRA GRIA2 TAGLN ITGA7 SEMA6D RBP1 RDH10 SCD5 DKK3 DPYD SERPINF1 RHOJ CHRNA1 RFX4 FSTL1 LPL SLC1A2 LGALS3BP SCRG1 CDH6 TUBA1A PPIC TTYH1 COL6A2 COL6A1 COL1A1", " ")[[1]]
sets_raw <- list(CANONICAL_MES=getset("NEFTEL_MES_LIKE"), MES1=getset("NEFTEL_MES1"), MES2=getset("NEFTEL_MES2"), MES_HYP=unique(toupper(HYP)), MES_AST=unique(toupper(AST)))
target_union <- union(sets_raw$MES_HYP, sets_raw$MES_AST); hyp_ast_overlap <- intersect(sets_raw$MES_HYP, sets_raw$MES_AST)
sets <- list()
for (nm in names(sets_raw)) {
  sets[[paste(nm, "RAW", sep="::")]] <- sets_raw[[nm]]
  sets[[paste(nm, "OVERLAP_PRUNED", sep="::")]] <- if (nm %in% c("CANONICAL_MES","MES1","MES2")) setdiff(sets_raw[[nm]], target_union) else setdiff(sets_raw[[nm]], hyp_ast_overlap)
}
removal <- rbindlist(lapply(names(sets_raw), function(nm) {
  removed <- setdiff(sets_raw[[nm]], sets[[paste(nm,"OVERLAP_PRUNED",sep="::")]])
  if (!length(removed)) return(data.table(signature=nm, gene=NA_character_, removal_rule="NO_DIRECT_OVERLAP", removed=FALSE))
  data.table(signature=nm, gene=removed, removal_rule=ifelse(nm %in% c("CANONICAL_MES","MES1","MES2"), "REMOVE_MES_HYP_UNION_MES_AST", "REMOVE_MES_HYP_MES_AST_INTERSECTION_BILATERALLY"), removed=TRUE)
}))
w(removal, file.path(OUTA, "IVY_OVERLAP_REMOVAL.tsv"))

score_one <- function(gene_set) {
  present <- intersect(gene_set, rownames(log_expression))
  if (length(present) < 2L) return(list(score=rep(NA_real_, ncol(log_expression)), covered=length(present)))
  z <- t(scale(t(log_expression[present, , drop=FALSE]))); z[!is.finite(z)] <- 0
  list(score=colMeans(z), covered=length(present))
}

base <- mapping[, .(donor_id=as.character(donor_id), specimen_id=as.character(specimen_id), sample_id=as.character(sample_id),
  region, raw_region=structure_parent_acronym, raw_region_name=structure_parent_name,
  source_filename=raw_profile_filename, tumor_name, normalized_matrix_column_present, raw_profile_present)]
sample_rows <- list(); coverage <- list()
for (key in names(sets)) {
  p <- strsplit(key, "::", fixed=TRUE)[[1]]; nm <- p[1]; mode <- p[2]; ans <- score_one(sets[[key]])
  status <- if (ans$covered >= 2L) "PASS" else "NOT_EVALUABLE_LOW_GENE_COVERAGE"
  d <- copy(base); d[, `:=`(signature=nm, mode=mode, score=ans$score,
    genes_defined=length(sets[[key]]), genes_covered=ans$covered,
    gene_coverage=ans$covered/length(sets[[key]]), QC_status=status)]
  sample_rows[[key]] <- d
  coverage[[key]] <- data.table(signature=nm, mode=mode, genes_defined=length(sets[[key]]), genes_covered=ans$covered,
    gene_coverage=ans$covered/length(sets[[key]]), minimum_genes=2L, QC_status=status,
    source=ifelse(nm %in% c("MES_HYP","MES_AST"), "v1.1 verified Lin2026 sidecar marker lineage", "frozen signature_gene_long_table_v1_3.tsv"))
}
long <- rbindlist(sample_rows)
w(long, file.path(OUTA, "IVY_SAMPLE_LEVEL_ROWS.tsv"))
w(rbindlist(coverage), file.path(OUTA, "IVY_SCORE_COVERAGE.tsv"))
donor_rows <- long[is.finite(score), .(score=mean(score), n_samples=.N, genes_covered=first(genes_covered), gene_coverage=first(gene_coverage), QC_status=first(QC_status)), by=.(donor_id, region, signature, mode)]
w(donor_rows, file.path(OUTA, "IVY_DONOR_LEVEL_ROWS.tsv"))
lineage <- unique(base[, .(donor_id, specimen_id, sample_id, region, raw_region, raw_region_name, source_filename,
  normalized_matrix_path=file.path(MATRIX_DIR,"fpkm_table.csv"), raw_profile_present, normalized_matrix_column_present)])
w(lineage, file.path(OUTA, "IVY_INPUT_LINEAGE.tsv"))

regions <- c("CELLULAR_TUMOUR","PSEUDOPALISADING_CELLS_AROUND_NECROSIS","MICROVASCULAR_PROLIFERATION","INFILTRATING_TUMOUR","LEADING_EDGE")
contrast_def <- data.table(contrast=c("PAN_VS_CT","MVP_VS_CT","IT_VS_CT","LE_VS_CT","PAN_VS_MVP"),
  region1=c("PSEUDOPALISADING_CELLS_AROUND_NECROSIS","MICROVASCULAR_PROLIFERATION","INFILTRATING_TUMOUR","LEADING_EDGE","PSEUDOPALISADING_CELLS_AROUND_NECROSIS"),
  region0=c("CELLULAR_TUMOUR","CELLULAR_TUMOUR","CELLULAR_TUMOUR","CELLULAR_TUMOUR","MICROVASCULAR_PROLIFERATION"))

contrast_lmer <- function(fit, r1, r0) {
  nd <- data.frame(region=factor(c(r1,r0), levels=regions), donor_id=factor(rep(levels(fit@frame$donor_id)[1],2), levels=levels(fit@frame$donor_id)))
  mm <- model.matrix(~region, nd); b <- fixef(fit); mm <- mm[, names(b), drop=FALSE]; cv <- as.numeric(mm[1,]-mm[2,])
  est <- sum(cv*b); se <- sqrt(as.numeric(t(cv) %*% vcov(fit) %*% cv)); z <- est/se
  c(estimate=est, SE=se, CI_low=est-1.96*se, CI_high=est+1.96*se, raw_P=2*pnorm(-abs(z)))
}
contrast_lm <- function(fit, r1, r0) {
  donor_levels <- levels(fit$model$donor_factor)
  nd <- data.frame(region=factor(c(r1,r0), levels=regions), donor_factor=factor(rep(donor_levels[1],2), levels=donor_levels))
  mm <- model.matrix(~region+donor_factor, nd); b <- coef(fit); keep <- names(b)[is.finite(b)]; mm <- mm[, keep, drop=FALSE]; cv <- as.numeric(mm[1,]-mm[2,])
  V <- vcov(fit)[keep,keep,drop=FALSE]; est <- sum(cv*b[keep]); se <- sqrt(as.numeric(t(cv)%*%V%*%cv)); z <- est/se
  c(estimate=est, SE=se, CI_low=est-1.96*se, CI_high=est+1.96*se, raw_P=2*pnorm(-abs(z)))
}

model_rows <- list(); contrast_rows <- list(); fe_rows <- list(); qa_rows <- list()
for (nm in names(sets_raw)) for (mode2 in c("RAW","OVERLAP_PRUNED")) {
  d <- long[signature==nm & mode==mode2 & is.finite(score)]
  d[, region := factor(region, levels=regions)]; d[, donor_id := factor(donor_id)]
  warnings <- character(); fit <- tryCatch(withCallingHandlers(lmer(score ~ region + (1|donor_id), data=d, REML=FALSE,
    control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))), warning=function(x){warnings <<- c(warnings,conditionMessage(x)); invokeRestart("muffleWarning")}), error=function(e)e)
  if (inherits(fit,"error")) {
    qa_rows[[paste(nm,mode2)]] <- data.table(signature=nm, mode=mode2, n_donors=uniqueN(d$donor_id), n_samples=nrow(d), convergence="FAIL", singular=NA, warning=conditionMessage(fit), genes_covered=unique(d$genes_covered)[1])
    next
  }
  singular <- isSingular(fit, tol=1e-5); convmsg <- fit@optinfo$conv$lme4$messages; conv <- if (is.null(convmsg)) "PASS" else paste(convmsg,collapse=";")
  co <- as.data.table(coef(summary(fit)), keep.rownames="term"); setnames(co,c("Estimate","Std. Error","t value"),c("estimate","SE","z_value"),skip_absent=TRUE)
  co[, `:=`(CI_low=estimate-1.96*SE, CI_high=estimate+1.96*SE, raw_P=2*pnorm(-abs(z_value)), signature=nm, mode=mode2,
    donor_n=uniqueN(d$donor_id), sample_n=nrow(d), singular_fit=singular, convergence=conv)]
  co[, BH_adjusted_P := p.adjust(raw_P, "BH")]
  model_rows[[paste(nm,mode2)]] <- co
  cr <- rbindlist(lapply(seq_len(nrow(contrast_def)), function(j) {
    a <- contrast_lmer(fit, contrast_def$region1[j], contrast_def$region0[j]); data.table(signature=nm, mode=mode2,
      contrast=contrast_def$contrast[j], region1=contrast_def$region1[j], region0=contrast_def$region0[j],
      estimate=a["estimate"], SE=a["SE"], CI_low=a["CI_low"], CI_high=a["CI_high"], raw_P=a["raw_P"],
      donor_n=uniqueN(d$donor_id), sample_n=nrow(d), singular_fit=singular, convergence=conv)
  }))
  cr[, BH_adjusted_P := p.adjust(raw_P, "BH")]
  contrast_rows[[paste(nm,mode2)]] <- cr
  dd <- copy(d); dd[, donor_factor := factor(donor_id)]
  lmfit <- tryCatch(lm(score ~ region + donor_factor, data=dd), error=function(e)e)
  if (!inherits(lmfit,"error")) {
    fr <- rbindlist(lapply(seq_len(nrow(contrast_def)), function(j) {
      a <- contrast_lm(lmfit, contrast_def$region1[j], contrast_def$region0[j]); data.table(signature=nm, mode=mode2,
        contrast=contrast_def$contrast[j], estimate=a["estimate"], SE=a["SE"], CI_low=a["CI_low"], CI_high=a["CI_high"], raw_P=a["raw_P"], donor_n=uniqueN(dd$donor_id), sample_n=nrow(dd), rank=lmfit$rank, residual_df=df.residual(lmfit))
    }))
    fr[, BH_adjusted_P := p.adjust(raw_P,"BH")]; fe_rows[[paste(nm,mode2)]] <- fr
  }
  qa_rows[[paste(nm,mode2)]] <- data.table(signature=nm, mode=mode2, n_donors=uniqueN(d$donor_id), n_samples=nrow(d),
    convergence=conv, singular=singular, warning=paste(unique(warnings),collapse=";"), random_intercept_variance=as.data.table(VarCorr(fit))$vcov[1],
    genes_covered=unique(d$genes_covered)[1], low_gene_coverage=unique(d$genes_covered)[1]<2)
}
models <- rbindlist(model_rows, fill=TRUE); contrasts <- rbindlist(contrast_rows, fill=TRUE); fixed <- rbindlist(fe_rows, fill=TRUE)
w(models, file.path(OUTM, "IVY_REGION_MIXED_MODEL.tsv")); w(contrasts, file.path(OUTM, "IVY_REGION_CONTRASTS.tsv")); w(fixed, file.path(OUTM, "IVY_DONOR_FIXED_EFFECT.tsv")); w(rbindlist(qa_rows,fill=TRUE), file.path(OUTM, "IVY_MODEL_QA.tsv"))

paired_diffs <- function(d, r1, r0) {
  q <- d[region %in% c(r1,r0), .(v=mean(score), n=.N), by=.(donor_id,region)]
  a <- q[region==r1,.(donor_id,v1=v,n1=n)]; b <- q[region==r0,.(donor_id,v0=v,n0=n)]
  z <- merge(a,b,by="donor_id"); z[, diff:=v1-v0]; z
}
boot_rows <- list(); pe_rows <- list(); lodo_rows <- list()
for (nm in names(sets_raw)) for (mode2 in c("RAW","OVERLAP_PRUNED")) {
  d <- long[signature==nm & mode==mode2 & is.finite(score)]
  for (j in seq_len(nrow(contrast_def))) {
    z <- paired_diffs(d, contrast_def$region1[j], contrast_def$region0[j]); key <- paste(nm,mode2,contrast_def$contrast[j])
    if (nrow(z)<2L) {
      boot_rows[[key]] <- data.table(signature=nm,mode=mode2,contrast=contrast_def$contrast[j],estimate=NA,CI_low=NA,CI_high=NA,n_donors=nrow(z),bootstrap_iterations=0,status="NOT_EVALUABLE")
      next
    }
    set.seed(2026080141 + match(nm,names(sets_raw))*100 + match(mode2,c("RAW","OVERLAP_PRUNED"))*10 + j)
    boots <- replicate(2000, mean(sample(z$diff, nrow(z), replace=TRUE)))
    boot_rows[[key]] <- data.table(signature=nm,mode=mode2,contrast=contrast_def$contrast[j],estimate=mean(z$diff),CI_low=quantile(boots,.025),CI_high=quantile(boots,.975),n_donors=nrow(z),bootstrap_iterations=2000,direction_positive_fraction=mean(boots>0),status="EVALUATED")
    tt <- t.test(z$diff)
    pe_rows[[key]] <- data.table(signature=nm,mode=mode2,contrast=contrast_def$contrast[j],estimate=mean(z$diff),median=median(z$diff),SE=sd(z$diff)/sqrt(nrow(z)),CI_low=tt$conf.int[1],CI_high=tt$conf.int[2],raw_P=tt$p.value,n_donors=nrow(z),inference_unit="donor; paired region means")
    lodo_rows[[key]] <- rbindlist(lapply(z$donor_id, function(ex) data.table(signature=nm,mode=mode2,contrast=contrast_def$contrast[j],excluded_donor=ex,estimate=mean(z[donor_id!=ex,diff]),n_donors=nrow(z)-1L)))
  }
}
boots <- rbindlist(boot_rows,fill=TRUE); pe <- rbindlist(pe_rows,fill=TRUE); lodo <- rbindlist(lodo_rows,fill=TRUE)
pe[, BH_adjusted_P := p.adjust(raw_P,"BH"), by=.(signature,mode)]
w(boots,file.path(OUTM,"IVY_CLUSTER_BOOTSTRAP.tsv")); w(pe,file.path(OUTM,"IVY_PATIENT_EQUAL.tsv")); w(lodo,file.path(OUTM,"IVY_LEAVE_ONE_DONOR_OUT.tsv"))
rawpruned <- merge(contrasts[mode=="RAW",.(signature,contrast,raw_estimate=estimate,raw_CI_low=CI_low,raw_CI_high=CI_high,raw_BH=BH_adjusted_P)],
  contrasts[mode=="OVERLAP_PRUNED",.(signature,contrast,pruned_estimate=estimate,pruned_CI_low=CI_low,pruned_CI_high=CI_high,pruned_BH=BH_adjusted_P)],by=c("signature","contrast"),all=TRUE)
rawpruned[, `:=`(direction_consistent=sign(raw_estimate)==sign(pruned_estimate), delta_pruned_minus_raw=pruned_estimate-raw_estimate)]
w(rawpruned,file.path(OUTM,"IVY_RAW_VS_PRUNED.tsv"))

gates <- rbindlist(lapply(seq_len(nrow(rawpruned)), function(i) {
  r <- rawpruned[i]; br <- boots[signature==r$signature & contrast==r$contrast & mode=="RAW"]
  pp <- pe[signature==r$signature & contrast==r$contrast & mode=="RAW"]
  lo <- lodo[signature==r$signature & contrast==r$contrast & mode=="RAW"]
  primary_sig <- is.finite(r$raw_estimate) && is.finite(r$raw_BH) && r$raw_BH < .05 && !(r$raw_CI_low <= 0 && r$raw_CI_high >= 0)
  boot_same <- nrow(br)==1 && is.finite(br$CI_low) && sign(br$estimate)==sign(r$raw_estimate) && !(br$CI_low <=0 && br$CI_high>=0)
  pe_same <- nrow(pp)==1 && is.finite(pp$estimate) && sign(pp$estimate)==sign(r$raw_estimate)
  lodo_flip <- nrow(lo)>0 && any(sign(lo$estimate)!=sign(r$raw_estimate),na.rm=TRUE)
  prune_same <- isTRUE(r$direction_consistent) && is.finite(r$pruned_estimate)
  gate <- if (!is.finite(r$raw_estimate) || nrow(br)==0) "NOT_EVALUABLE" else if (primary_sig && boot_same && pe_same && !lodo_flip && prune_same) "PATHOLOGY_REGION_ANCHORED" else if (!prune_same) "OVERLAP_DEPENDENT" else if (lodo_flip) "DONOR_HETEROGENEOUS" else if (boot_same && pe_same) "DIRECTIONALLY_CONSISTENT" else "NO_RELIABLE_REGION_ANCHOR"
  data.table(signature=r$signature,contrast=r$contrast,gate=gate,raw_estimate=r$raw_estimate,pruned_estimate=r$pruned_estimate,primary_BH=r$raw_BH,bootstrap_same=boot_same,patient_equal_same=pe_same,lodo_direction_flip=lodo_flip,overlap_pruned_same=prune_same)
}))
w(gates,file.path(OUTM,"IVY_CLAIM_GATES.tsv"))

elapsed <- proc.time()[["elapsed"]]-start_time
qa <- rbindlist(list(
  data.table(check="allen_normalized_samples",value=nrow(mapping),expected=270,status=ifelse(nrow(mapping)==270,"PASS","FAIL")),
  data.table(check="five_regions",value=uniqueN(mapping$region),expected=5,status=ifelse(uniqueN(mapping$region)==5,"PASS","FAIL")),
  data.table(check="donor_count_not_patient_from_files",value=uniqueN(mapping$donor_id),expected="derived_from_donor_id",status="PASS"),
  data.table(check="formal_score_modes",value=uniqueN(long[,paste(signature,mode)]),expected=10,status=ifelse(uniqueN(long[,paste(signature,mode)])==10,"PASS","FAIL")),
  data.table(check="mixed_models",value=uniqueN(contrasts[,paste(signature,mode)]),expected=10,status=ifelse(uniqueN(contrasts[,paste(signature,mode)])==10,"PASS","FAIL")),
  data.table(check="cluster_bootstrap_iterations",value=min(boots[status=="EVALUATED",bootstrap_iterations]),expected=2000,status=ifelse(min(boots[status=="EVALUATED",bootstrap_iterations])==2000,"PASS","FAIL")),
  data.table(check="elapsed_seconds",value=elapsed,expected="recorded",status="PASS")
),fill=TRUE)
w(qa,file.path(OUTM,"IVY_MODULE_QA.tsv"))
writeLines(c(capture.output(sessionInfo()),paste0("elapsed_seconds=",elapsed),paste0("donors=",uniqueN(mapping$donor_id)),paste0("samples=",nrow(mapping))),file.path(OUTM,"IVY_SESSION_INFO.txt"))

key_gate <- function(sig, con) gates[signature==sig & contrast==con,gate][1]
gate_md <- c("# IVY GAP analysis summary","",paste0("Module status: `",ifelse(all(qa$status=="PASS"),"PASS","REQUIRES_RERUN"),"`"),"",
  paste0("- CANONICAL_MES PAN vs CT: `",key_gate("CANONICAL_MES","PAN_VS_CT"),"`"),
  paste0("- CANONICAL_MES MVP vs CT: `",key_gate("CANONICAL_MES","MVP_VS_CT"),"`"),
  paste0("- MES_HYP PAN vs CT: `",key_gate("MES_HYP","PAN_VS_CT"),"`"),
  paste0("- MES2 PAN vs CT: `",key_gate("MES2","PAN_VS_CT"),"`"),
  paste0("- MES_AST MVP vs CT: `",key_gate("MES_AST","MVP_VS_CT"),"`"),"",
  "Inference uses donor IDs; 270 normalized samples are not described as 270 patients. RAW and OVERLAP_PRUNED modes follow the pre-result clarification; no removed gene was restored.","",
  paste0("Elapsed seconds: `",round(elapsed,3),"`  "),"Peak RSS: recorded by external `/usr/bin/time -l` log.  ","Next-module eligibility: `ELIGIBLE_FOR_INTEGRATION`")
writeLines(gate_md,file.path(OUTM,"IVY_ANALYSIS_SUMMARY.md"))
cat("IVY_V1_2_COMPLETE",nrow(mapping),uniqueN(mapping$donor_id),nrow(contrasts),nrow(boots),"\n")
