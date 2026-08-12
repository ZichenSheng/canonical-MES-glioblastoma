#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
suppressPackageStartupMessages({library(data.table); library(singscore); library(metafor)})

DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset="data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset="results")
REPO_ROOT <- Sys.getenv("GBM_MES_REPO_ROOT", unset=".")
RUN <- file.path(RESULT_ROOT,"interpretability")
OUT <- file.path(RUN,"state_specificity")
V1 <- file.path(DATA_ROOT,"prepared","core")
V11 <- file.path(RESULT_ROOT,"scoring")
SIG <- file.path(REPO_ROOT,"resources","signatures","signature_gene_sets.tsv")
set.seed(as.integer(20260801101 %% .Machine$integer.max))

w <- function(x, n) fwrite(x, file.path(OUT, n), sep="\t", quote=FALSE, na="NA")
zvec <- function(x) { s <- sd(x, na.rm=TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x))); (x-mean(x,na.rm=TRUE))/s }
score_zmean <- function(expr, genes) {
  g <- intersect(genes, rownames(expr)); if (length(g) < 2) return(rep(NA_real_, ncol(expr)))
  x <- expr[g,,drop=FALSE]; mu <- rowMeans(x); ss <- apply(x,1,sd); ok <- is.finite(ss) & ss > 0
  if (sum(ok) < 2) return(rep(NA_real_, ncol(expr)))
  colMeans((x[ok,,drop=FALSE]-mu[ok])/ss[ok])
}
context_score <- function(d) rowMeans(as.data.frame(lapply(as.data.frame(d)[,c("hypoxia","EMT","matrix","myeloid"),drop=FALSE], zvec)), na.rm=TRUE)
boot_pair <- function(target, control, context, B=5000) {
  ok <- complete.cases(target,control,context); target <- target[ok]; control <- control[ok]; context <- context[ok]; n <- length(target)
  observed_target <- suppressWarnings(cor(target,context,method="spearman")); observed_control <- suppressWarnings(cor(control,context,method="spearman"))
  b_target <- b_control <- b_delta <- rep(NA_real_,B)
  for (b in seq_len(B)) { ii <- sample.int(n,n,TRUE); b_target[b] <- suppressWarnings(cor(target[ii],context[ii],method="spearman")); b_control[b] <- suppressWarnings(cor(control[ii],context[ii],method="spearman")); b_delta[b] <- b_target[b]-b_control[b] }
  list(n=n, rho_target=observed_target, rho_control=observed_control, delta=observed_target-observed_control,
       target_ci=quantile(b_target,c(.025,.975),na.rm=TRUE,names=FALSE), control_ci=quantile(b_control,c(.025,.975),na.rm=TRUE,names=FALSE),
       delta_ci=quantile(b_delta,c(.025,.975),na.rm=TRUE,names=FALSE), delta_se=sd(b_delta,na.rm=TRUE), delta_p=2*min(mean(b_delta<=0,na.rm=TRUE),mean(b_delta>=0,na.rm=TRUE)))
}
cor_row <- function(cohort, state, method, x, y, B=5000) {
  ok <- complete.cases(x,y); x <- x[ok]; y <- y[ok]; n <- length(x); rho <- suppressWarnings(cor(x,y,method="spearman"));
  bs <- replicate(B,{ii<-sample.int(n,n,TRUE);suppressWarnings(cor(x[ii],y[ii],method="spearman"))})
  ct <- suppressWarnings(cor.test(x,y,method="spearman",exact=FALSE)); ci <- quantile(bs,c(.025,.975),na.rm=TRUE)
  data.table(cohort,state,method,n,rho,p_value=ct$p.value,ci_low=ci[1],ci_high=ci[2],genes_retained=NA_integer_)
}
meta_cor <- function(d, endpoint) {
  keep <- is.finite(d$rho)&d$n>3; d <- d[keep]; z <- atanh(pmax(pmin(d$rho,.999999),-.999999)); vi <- 1/(d$n-3)
  fit <- rma.uni(z,vi=vi,method="REML"); pr <- predict(fit)
  data.table(endpoint,k=nrow(d),estimate=tanh(as.numeric(fit$b)),ci_low=tanh(fit$ci.lb),ci_high=tanh(fit$ci.ub),prediction_low=tanh(pr$pi.lb),prediction_high=tanh(pr$pi.ub),I2=fit$I2,tau2=fit$tau2,direction_consistency=paste(sign(d$rho),collapse=";"))
}
meta_delta <- function(d, control) {
  d <- d[is.finite(delta)&is.finite(delta_se)&delta_se>0]; fit <- rma.uni(yi=d$delta,vi=d$delta_se^2,method="REML"); pr <- predict(fit)
  data.table(control,k=nrow(d),delta_rho=as.numeric(fit$b),ci_low=fit$ci.lb,ci_high=fit$ci.ub,prediction_low=pr$pi.lb,prediction_high=pr$pi.ub,I2=fit$I2,tau2=fit$tau2,direction_consistency=paste(sign(d$delta),collapse=";"))
}

reg <- fread(SIG); reg[,gene:=toupper(trimws(gene_symbol_clean))]
map <- c(CANONICAL_MES="NEFTEL_MES_LIKE", AC_LIKE="NEFTEL_AC_LIKE", OPC_LIKE="NEFTEL_OPC_LIKE", NPC_LIKE="NEFTEL_NPC_LIKE", PROLIFERATION="PROLIFERATION_CELL_CYCLE", hypoxia="HALLMARK_HYPOXIA", EMT="HALLMARK_EMT", matrix="NABA_CORE_MATRISOME", myeloid="MYELOID_CORE_IDENTITY")
sets <- lapply(map,function(id) unique(reg[signature_id==id & !is.na(gene) & nzchar(gene),gene]))
controls <- c("AC_LIKE","OPC_LIKE","NPC_LIKE","PROLIFERATION")

cohort_members <- fread(file.path(V1,"03_cohort/PATIENT_LEVEL_COHORT.tsv"))
tcga_members <- fread(file.path(V11,"05_tcga_strict/TCGA_STRICT_COHORT.tsv"))
paths <- list(
  CGGA325=file.path(DATA_ROOT,"bulk","cgga325_expression.tsv"),
  CGGA693=file.path(DATA_ROOT,"bulk","cgga693_expression.tsv"),
  TCGA=file.path(DATA_ROOT,"bulk","tcga_gbm_expression.tsv.gz")
)
raw_results <- list(); pruned_results <- list(); sing_results <- list(); deltas <- list()

for (cohort_name in names(paths)) {
  if (cohort_name=="TCGA") { ids <- tcga_members[strict_inclusion==TRUE,sample_id] } else { ids <- cohort_members[cohort==cohort_name & strict_included==TRUE,patient_id] }
  hdr <- names(fread(paths[[cohort_name]],nrows=0,check.names=FALSE)); stopifnot(all(ids %in% hdr))
  dt <- fread(paths[[cohort_name]],select=c(hdr[1],ids),check.names=FALSE); genes <- toupper(trimws(dt[[1]])); expr <- as.matrix(dt[,-1]); storage.mode(expr)<-"double"; rownames(expr)<-genes
  if (anyDuplicated(genes)) { u<-unique(genes); expr<-rowsum(expr,genes,reorder=FALSE)/as.vector(table(factor(genes,levels=u))) }
  colnames(expr)<-ids; if (cohort_name!="TCGA") expr <- log2(expr+1)
  sc <- data.table(patient_id=ids)
  for (nm in names(sets)) sc[[nm]] <- score_zmean(expr,sets[[nm]])
  sc[,context:=context_score(sc)]
  for (state in c("CANONICAL_MES",controls)) raw_results[[paste(cohort_name,state)]] <- cor_row(cohort_name,state,"mean_z",sc[[state]],sc$context)

  rd <- rankGenes(expr); sing <- data.table(patient_id=ids)
  for (nm in names(sets)) sing[[nm]] <- simpleScore(rd,upSet=intersect(sets[[nm]],rownames(expr)),centerScore=TRUE,knownDirection=TRUE)$TotalScore
  sing[,context:=context_score(sing)]
  for (state in c("CANONICAL_MES",controls)) sing_results[[paste(cohort_name,state)]] <- cor_row(cohort_name,state,"singscore",sing[[state]],sing$context)

  ctx_union <- unique(unlist(sets[c("hypoxia","EMT","matrix","myeloid")]))
  for (state in c("CANONICAL_MES",controls)) {
    state_genes <- setdiff(sets[[state]],ctx_union); ctx_sets <- lapply(sets[c("hypoxia","EMT","matrix","myeloid")],setdiff,y=sets[[state]])
    ps <- data.table(state=score_zmean(expr,state_genes)); for (nm in names(ctx_sets)) ps[[nm]] <- score_zmean(expr,ctx_sets[[nm]]); ps[,context:=context_score(ps)]
    rr <- cor_row(cohort_name,state,"mean_z_bilateral_pruned",ps$state,ps$context); rr[,genes_retained:=length(intersect(state_genes,rownames(expr)))]; pruned_results[[paste(cohort_name,state)]] <- rr
  }
  for (control in controls) {
    bp <- boot_pair(sc$CANONICAL_MES,sc[[control]],sc$context)
    deltas[[paste(cohort_name,control)]] <- data.table(cohort=cohort_name,control,n=bp$n,rho_canonical=bp$rho_target,rho_control=bp$rho_control,delta=bp$delta,ci_low=bp$delta_ci[1],ci_high=bp$delta_ci[2],delta_se=bp$delta_se,p_value=bp$delta_p,bootstrap_replicates=5000L,paired_patients=TRUE)
  }
}

raw <- rbindlist(raw_results); pruned <- rbindlist(pruned_results); sing <- rbindlist(sing_results); delta <- rbindlist(deltas)
raw[,FDR:=p.adjust(p_value,"BH")]; pruned[,FDR:=p.adjust(p_value,"BH")]; sing[,FDR:=p.adjust(p_value,"BH")]; delta[,FDR:=p.adjust(p_value,"BH")]
meta <- rbindlist(lapply(split(raw,raw$state),function(x)meta_cor(x,paste0(unique(x$state),"_mean_z"))))
delta_meta <- rbindlist(lapply(controls,function(x)meta_delta(delta[control==x],x)))
w(raw,"STRICT_BULK_CONTROL_RESULTS.tsv")
w(rbind(meta,delta_meta,fill=TRUE),"STRICT_BULK_CONTROL_META.tsv")
w(rbind(delta,delta_meta,fill=TRUE),"STRICT_BULK_DEPENDENT_DELTA.tsv")
w(rbind(cbind(analysis="RAW",raw),cbind(analysis="BILATERAL_OVERLAP_PRUNED",pruned),fill=TRUE),"STRICT_BULK_RAW_VS_PRUNED.tsv")
w(rbind(cbind(analysis="MEAN_Z",raw),cbind(analysis="SINGSCORE",sing),fill=TRUE),"STRICT_BULK_SCORING_SENSITIVITY.tsv")

bench <- fread(file.path(V11,"04_pseudobulk_benchmark/SCORE_INTERPRETABILITY_BENCHMARK.tsv"))
bench <- bench[score_name %chin% c("AC_like","OPC_like","NPC_like","proliferation")]
bench[,state:=fcase(score_name=="AC_like","AC_LIKE",score_name=="OPC_like","OPC_LIKE",score_name=="NPC_like","NPC_LIKE",default="PROLIFERATION")]
common_cols <- c("experiment","compartment","state","composition_slope","malignant_state_slope","interaction","patient_level_variability","study_level_variability","model_R2","patients","studies","pseudobulks","formal_result_role")
w(bench[experiment=="A_FIXED_MALIGNANT",..common_cols],"PSEUDOBULK_CONTROL_COMPOSITION.tsv")
w(bench[experiment=="B_FIXED_COMPOSITION",..common_cols],"PSEUDOBULK_CONTROL_OFFTARGET.tsv")
w(bench[experiment=="C_TWO_DIMENSIONAL",..common_cols],"PSEUDOBULK_CONTROL_INTERACTION.tsv")
writeLines(capture.output(sessionInfo()),file.path(RUN,"session_info/MODULE1_BULK_PSEUDOBULK_SESSION_INFO.txt"))
