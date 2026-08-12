#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
suppressPackageStartupMessages({library(data.table); library(MASS); library(splines)})

DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset="data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset="results")
REPO_ROOT <- Sys.getenv("GBM_MES_REPO_ROOT", unset=".")
RUN <- file.path(RESULT_ROOT,"patient_realization")
V10 <- file.path(DATA_ROOT,"prepared","core")
V11 <- file.path(RESULT_ROOT,"scoring")
V13 <- file.path(RESULT_ROOT,"interpretability")
V14 <- file.path(RESULT_ROOT,"provenance")
SIG <- file.path(REPO_ROOT,"resources","signatures","signature_gene_sets.tsv")
EXPR_PATHS <- c(
  CGGA325=file.path(DATA_ROOT,"bulk","cgga325_expression.tsv"),
  CGGA693=file.path(DATA_ROOT,"bulk","cgga693_expression.tsv")
)
SEEDS <- list(axis_boot=2026080911L, axis_perm=2026080912L, iso_boot=2026080921L,
              nn_perm=2026080922L, state_boot=2026080931L, state_perm=2026080932L)
B <- 5000L
THRESHOLDS <- c(0.05,0.10,0.15,0.20)

write_tab <- function(x, rel) {
  path <- file.path(RUN,rel); dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  fwrite(x,path,sep="\t",quote=FALSE,na="NA")
}
zv <- function(x) {
  s <- sd(x,na.rm=TRUE)
  if(!is.finite(s)||s==0) return(rep(NA_real_,length(x)))
  (x-mean(x,na.rm=TRUE))/s
}
qtiles <- function(x) unname(quantile(x,c(.25,.5,.75,.9),na.rm=TRUE,type=7))
score_zmean <- function(expr,genes) {
  present <- intersect(genes,rownames(expr))
  if(length(present)<2L) return(list(score=rep(NA_real_,ncol(expr)),covered=length(present),defined=length(genes)))
  x <- expr[present,,drop=FALSE]; m <- rowMeans(x); s <- apply(x,1,sd)
  ok <- is.finite(s)&s>0
  if(sum(ok)<2L) return(list(score=rep(NA_real_,ncol(expr)),covered=sum(ok),defined=length(genes)))
  list(score=colMeans((x[ok,,drop=FALSE]-m[ok])/s[ok]),covered=sum(ok),defined=length(genes))
}
sha <- function(path) unname(tools::md5sum(path)) # display-only fallback; true SHA sources come from governance inventory

cat("FORMAL_ANALYSIS_START",format(Sys.time(),tz="UTC",usetz=TRUE),"\n")

# Frozen memberships and inputs.
mem <- rbindlist(lapply(c("CGGA325","CGGA693"),function(co) {
  fread(file.path(RUN,"00_governance/freeze",paste0("membership_",tolower(co),"_patients_v1_6.tsv")))
}))
stopifnot(mem[cohort=="CGGA325",uniqueN(patient_id)]==73L,mem[cohort=="CGGA693",uniqueN(patient_id)]==108L)
stopifnot(!anyDuplicated(mem[,.(cohort,patient_id)]))

bulk <- fread(file.path(V10,"05_bulk/BULK_PATIENT_SCORES.tsv"))
bpstate <- fread(file.path(V13,"01_state_controls/BAYESPRISM_STATE_PATIENT_SCORES_INTERNAL.tsv"))
frac <- fread(file.path(V11,"07_bayesprism/CELL_FRACTION_ESTIMATES_REUSED.tsv"))
malmes <- fread(file.path(V11,"07_bayesprism/MALIGNANT_EXPRESSION_ESTIMATES.tsv"))
reg <- fread(SIG); reg[,gene:=toupper(trimws(gene_symbol_clean))]
reg <- reg[!is.na(gene)&nzchar(gene)]
setmap <- c(CANONICAL_MES="NEFTEL_MES_LIKE",AC_LIKE="NEFTEL_AC_LIKE",OPC_LIKE="NEFTEL_OPC_LIKE",
            NPC_LIKE="NEFTEL_NPC_LIKE",PROLIFERATION="PROLIFERATION_CELL_CYCLE",
            HYPOXIA="HALLMARK_HYPOXIA",MATRIX="NABA_CORE_MATRISOME")
sets <- lapply(setmap,function(s)unique(reg[signature_id==s,gene]))
expected_sizes <- c(CANONICAL_MES=95L,AC_LIKE=39L,OPC_LIKE=50L,NPC_LIKE=89L,PROLIFERATION=327L,HYPOXIA=200L,MATRIX=275L)
stopifnot(identical(as.integer(lengths(sets)),as.integer(expected_sizes[names(sets)])))

# Frozen-v1.3 compatible ecology rescoring, including state-specific bilateral pruning.
eco_list <- list(); eco_cov <- list(); expr_official_qa <- list()
state_names <- c("CANONICAL_MES","AC_LIKE","OPC_LIKE","NPC_LIKE","PROLIFERATION")
for(ci in seq_along(EXPR_PATHS)) {
  co <- names(EXPR_PATHS)[ci]; ids <- mem[cohort==co,patient_id]
  hdr <- names(fread(EXPR_PATHS[[ci]],nrows=0)); miss <- setdiff(ids,hdr)
  if(length(miss)) stop("Expression membership mismatch: ",co," ",paste(miss,collapse=","))
  dt <- fread(EXPR_PATHS[[ci]],select=c(hdr[1],ids),check.names=FALSE)
  genes <- toupper(trimws(dt[[1]])); mat <- as.matrix(dt[,-1]); storage.mode(mat)<-"double"; rownames(mat)<-genes
  if(anyDuplicated(genes)) {
    u <- unique(genes); mat <- rowsum(mat,genes,reorder=FALSE)/as.vector(table(factor(genes,levels=u)))
  }
  expr <- log2(mat+1); colnames(expr)<-ids
  e <- data.table(cohort=co,patient_id=ids)
  raw_h <- score_zmean(expr,sets$HYPOXIA); raw_x <- score_zmean(expr,sets$MATRIX)
  e[,`:=`(hypoxia_raw=raw_h$score,matrix_raw=raw_x$score)]
  eco_cov[[length(eco_cov)+1L]] <- data.table(cohort=co,state="RAW",axis=c("HYPOXIA","MATRIX"),genes_defined=c(raw_h$defined,raw_x$defined),genes_covered=c(raw_h$covered,raw_x$covered),genes_removed_direct_overlap=0L)
  for(st in state_names) {
    hp <- setdiff(sets$HYPOXIA,sets[[st]]); xp <- setdiff(sets$MATRIX,sets[[st]])
    hs <- score_zmean(expr,hp); xs <- score_zmean(expr,xp)
    e[[paste0("hypoxia_pruned_",st)]] <- hs$score
    e[[paste0("matrix_pruned_",st)]] <- xs$score
    eco_cov[[length(eco_cov)+1L]] <- data.table(cohort=co,state=st,axis=c("HYPOXIA","MATRIX"),genes_defined=c(length(hp),length(xp)),genes_covered=c(hs$covered,xs$covered),genes_removed_direct_overlap=c(length(intersect(sets$HYPOXIA,sets[[st]])),length(intersect(sets$MATRIX,sets[[st]]))))
  }
  recomputed <- score_zmean(expr,sets$CANONICAL_MES)$score
  official <- bulk[cohort==co][match(ids,patient_id),canonical_MES]
  expr_official_qa[[length(expr_official_qa)+1L]] <- data.table(cohort=co,n=length(ids),max_abs_difference=max(abs(recomputed-official)),MAE=mean(abs(recomputed-official)),pearson=cor(recomputed,official),spearman=cor(recomputed,official,method="spearman"),status=ifelse(max(abs(recomputed-official))<1e-12,"PASS_EXACT","FAIL"))
  eco_list[[co]] <- e
}
eco <- rbindlist(eco_list)
write_tab(rbindlist(eco_cov),"01_reference/ECOLOGY_PRUNING_GENE_COVERAGE.tsv")
write_tab(rbindlist(expr_official_qa),"01_reference/OFFICIAL_SCORE_RECOMPUTATION_QA.tsv")
stopifnot(all(rbindlist(expr_official_qa)$status=="PASS_EXACT"))

# Master patient reference.
d <- merge(mem[,.(cohort,patient_id)],bulk[,.(cohort,patient_id,bulk_mes_raw=canonical_MES)],by=c("cohort","patient_id"),all.x=TRUE)
d <- merge(d,bpstate[,.(cohort,patient_id,malignant_mes_raw=malignant_CANONICAL_MES,myeloid_fraction,vascular_stromal_fraction,
                        malignant_AC_LIKE,malignant_OPC_LIKE,malignant_NPC_LIKE,malignant_PROLIFERATION,
                        coverage_CANONICAL_MES,coverage_fraction_CANONICAL_MES,coverage_AC_LIKE,coverage_fraction_AC_LIKE,
                        coverage_OPC_LIKE,coverage_fraction_OPC_LIKE,coverage_NPC_LIKE,coverage_fraction_NPC_LIKE,
                        coverage_PROLIFERATION,coverage_fraction_PROLIFERATION,
                        bulk_AC_LIKE,bulk_OPC_LIKE,bulk_NPC_LIKE,bulk_PROLIFERATION,posterior_median_cv)],by=c("cohort","patient_id"),all.x=TRUE)
d <- merge(d,eco,by=c("cohort","patient_id"),all.x=TRUE)
stopifnot(nrow(d)==181L,!anyDuplicated(d[,.(cohort,patient_id)]))
d[,`:=`(hypoxia_pruned=hypoxia_pruned_CANONICAL_MES,matrix_pruned=matrix_pruned_CANONICAL_MES)]
d[,`:=`(bulk_mes_z=zv(bulk_mes_raw),malignant_mes_z=zv(malignant_mes_raw),myeloid_z=zv(myeloid_fraction),
        hypoxia_z=zv(hypoxia_pruned),matrix_z=zv(matrix_pruned),vascular_z=zv(vascular_stromal_fraction),
        hypoxia_raw_z=zv(hypoxia_raw),matrix_raw_z=zv(matrix_raw)),by=cohort]
q4_rows <- list()
for(co in c("CGGA325","CGGA693")) {
  threshold <- quantile(d[cohort==co,bulk_mes_raw],.75,type=7)
  d[cohort==co,`:=`(mes_q4=bulk_mes_raw>=threshold,mes_q4_threshold=threshold)]
  q4_rows[[co]] <- data.table(cohort=co,n=.N,q4_threshold=threshold,q4_n=d[cohort==co,sum(mes_q4)])
}
q4 <- rbindlist(q4_rows)
if(!identical(q4$q4_n,c(19L,27L))) stop("Q4 count mismatch: ",paste(q4$q4_n,collapse=","))

inv <- fread(file.path(RUN,"00_governance/freeze/source_inventory.tsv"))
src <- setNames(inv$path,inv$source_id); shv <- setNames(inv$sha256,inv$source_id)
d[,`:=`(
  hypoxia_gene_coverage=rbindlist(eco_cov)[state=="CANONICAL_MES"&axis=="HYPOXIA",genes_covered][match(co,rbindlist(eco_cov)[state=="CANONICAL_MES"&axis=="HYPOXIA",cohort])],
  matrix_gene_coverage=rbindlist(eco_cov)[state=="CANONICAL_MES"&axis=="MATRIX",genes_covered][match(co,rbindlist(eco_cov)[state=="CANONICAL_MES"&axis=="MATRIX",cohort])]
),by=.(co=cohort)]
d[,missing_reason:=ifelse(rowSums(!is.finite(as.matrix(.SD)))>0,"PRIMARY_AXIS_MISSING","NONE"),.SDcols=c("bulk_mes_raw","malignant_mes_raw","myeloid_fraction","hypoxia_pruned","matrix_pruned","vascular_stromal_fraction")]
d[,coverage_flags:=paste0("MALIGNANT_CANONICAL=",coverage_CANONICAL_MES,"/95;HYPOXIA_PRUNED=",hypoxia_gene_coverage,"/",rbindlist(eco_cov)[state=="CANONICAL_MES"&axis=="HYPOXIA",genes_defined][match(cohort,rbindlist(eco_cov)[state=="CANONICAL_MES"&axis=="HYPOXIA",cohort])],";MATRIX_PRUNED=",matrix_gene_coverage,"/",rbindlist(eco_cov)[state=="CANONICAL_MES"&axis=="MATRIX",genes_defined][match(cohort,rbindlist(eco_cov)[state=="CANONICAL_MES"&axis=="MATRIX",cohort])])]
d[,all_source_paths:=paste(c(src[c("SRC-V10-BULK-SCORES","SRC-V13-STATE-SCORES","SRC-V11-CELL-FRACTIONS","SRC-CGGA325-EXPRESSION","SRC-CGGA693-EXPRESSION","SRC-SIGNATURE-REGISTRY")]),collapse=";")]
d[,all_source_sha256:=paste(c(shv[c("SRC-V10-BULK-SCORES","SRC-V13-STATE-SCORES","SRC-V11-CELL-FRACTIONS","SRC-CGGA325-EXPRESSION","SRC-CGGA693-EXPRESSION","SRC-SIGNATURE-REGISTRY")]),collapse=";")]
d[,`:=`(bulk_source_path=src[["SRC-V10-BULK-SCORES"]],bulk_source_sha256=shv[["SRC-V10-BULK-SCORES"]],
        malignant_source_path=src[["SRC-V13-STATE-SCORES"]],malignant_source_sha256=shv[["SRC-V13-STATE-SCORES"]],
        fraction_source_path=src[["SRC-V11-CELL-FRACTIONS"]],fraction_source_sha256=shv[["SRC-V11-CELL-FRACTIONS"]],
        ecology_source_path=ifelse(cohort=="CGGA325",src[["SRC-CGGA325-EXPRESSION"]],src[["SRC-CGGA693-EXPRESSION"]]),
        ecology_source_sha256=ifelse(cohort=="CGGA325",shv[["SRC-CGGA325-EXPRESSION"]],shv[["SRC-CGGA693-EXPRESSION"]]),
        signature_source_path=src[["SRC-SIGNATURE-REGISTRY"]],signature_source_sha256=shv[["SRC-SIGNATURE-REGISTRY"]])]
if(any(d$missing_reason!="NONE")) stop("Primary axis missing")
master_cols <- c("cohort","patient_id","bulk_mes_raw","bulk_mes_z","malignant_mes_raw","malignant_mes_z","myeloid_fraction","myeloid_z",
                 "hypoxia_pruned","hypoxia_z","matrix_pruned","matrix_z","vascular_stromal_fraction","vascular_z","hypoxia_raw","hypoxia_raw_z","matrix_raw","matrix_raw_z",
                 "mes_q4","mes_q4_threshold","coverage_CANONICAL_MES","hypoxia_gene_coverage","matrix_gene_coverage","missing_reason","coverage_flags",
                 "bulk_source_path","bulk_source_sha256","malignant_source_path","malignant_source_sha256","fraction_source_path","fraction_source_sha256",
                 "ecology_source_path","ecology_source_sha256","signature_source_path","signature_source_sha256","all_source_paths","all_source_sha256")
setorder(d,cohort,patient_id); write_tab(d[,.SD,.SDcols=master_cols],"01_reference/PATIENT_REALIZATION_MASTER.tsv")
write_tab(q4,"01_reference/Q4_DENOMINATOR_QA.tsv")

# Prediction functions.
linear_lopo <- function(x,y) {
  fit <- lm(y~x); h <- hatvalues(fit); e <- residuals(fit); den <- 1-h
  pred <- y-e/den; pred[!is.finite(pred)] <- NA_real_; pred
}
linear_lopo_cluster <- function(x,y,key=seq_along(y)) {
  n<-length(y); pred<-rep(NA_real_,n); ux<-unique(key)
  N<-n; sx<-sum(x); sy<-sum(y); sxx<-sum(x*x); sxy<-sum(x*y)
  for(k in ux){ii<-which(key==k);w<-length(ii);nt<-N-w;sxt<-sx-sum(x[ii]);syt<-sy-sum(y[ii]);sxxt<-sxx-sum(x[ii]^2);sxyt<-sxy-sum(x[ii]*y[ii]);den<-sxxt-sxt*sxt/nt
    if(is.finite(den)&&abs(den)>.Machine$double.eps){beta<-(sxyt-sxt*syt/nt)/den;alpha<-(syt-beta*sxt)/nt;pred[ii]<-alpha+beta*x[ii]}}
  pred
}
model_lopo <- function(x,y,model) {
  if(model=="LINEAR") return(linear_lopo(x,y))
  out <- rep(NA_real_,length(y))
  for(i in seq_along(y)) {
    tr <- data.frame(y=y[-i],x=x[-i]); nd <- data.frame(x=x[i])
    out[i] <- tryCatch({
      if(model=="HUBER") predict(MASS::rlm(y~x,data=tr,maxit=100),newdata=nd)
      else predict(lm(y~splines::ns(x,df=3),data=tr),newdata=nd)
    },error=function(e) NA_real_)
  }
  out
}
pred_metrics <- function(x,y,pred,model) {
  ok <- is.finite(x)&is.finite(y)&is.finite(pred); xx<-x[ok]; yy<-y[ok]; pp<-pred[ok]
  den <- sum((yy-mean(yy))^2)
  ins <- if(model=="LINEAR") summary(lm(yy~xx))$r.squared else if(model=="HUBER") cor(yy,predict(MASS::rlm(yy~xx)))^2 else summary(lm(yy~splines::ns(xx,df=3)))$r.squared
  c(pearson=cor(xx,yy),spearman=cor(xx,yy,method="spearman"),in_sample_R2=ins,Q2=1-sum((yy-pp)^2)/den,MAE=mean(abs(yy-pp)),RMSE=sqrt(mean((yy-pp)^2)))
}
axis_raw_cols <- c(MALIGNANT_STATE="malignant_mes_raw",MYELOID="myeloid_fraction",HYPOXIA="hypoxia_pruned",MATRIX="matrix_pruned",VASCULAR_STROMAL="vascular_stromal_fraction")
axis_z_cols <- c(MALIGNANT_STATE="malignant_mes_z",MYELOID="myeloid_z",HYPOXIA="hypoxia_z",MATRIX="matrix_z",VASCULAR_STROMAL="vascular_z")
pred_rows <- list(); conditioned <- list(); pred_boot_internal <- list(); ri<-0L; bi<-0L
for(coi in seq_along(c("CGGA325","CGGA693"))) {
  co <- c("CGGA325","CGGA693")[coi]; xdat <- copy(d[cohort==co]); x <- xdat$bulk_mes_z
  rdt <- xdat[,.(cohort,patient_id,bulk_mes_z)]
  for(ai in seq_along(axis_z_cols)) {
    axis <- names(axis_z_cols)[ai]; y <- xdat[[axis_z_cols[[ai]]]]
    linpred <- linear_lopo(x,y); resid <- y-linpred; rz <- zv(resid)
    rdt[[paste0("residual_raw_",axis)]] <- resid; rdt[[paste0("residual_z_",axis)]] <- rz; rdt[[paste0("lopo_prediction_",axis)]] <- linpred
    set.seed(SEEDS$axis_perm + coi*100L + ai)
    obs_p <- cor(x,y); obs_s <- cor(x,y,method="spearman")
    perm_p <- replicate(B,cor(x,sample(y),method="pearson")); perm_s <- replicate(B,cor(x,sample(y),method="spearman"))
    p_p <- (1+sum(abs(perm_p)>=abs(obs_p)))/(B+1); p_s <- (1+sum(abs(perm_s)>=abs(obs_s)))/(B+1)
    set.seed(SEEDS$axis_boot + coi*100L + ai)
    qb <- replicate(B,{ii<-sample.int(length(y),length(y),TRUE); pp<-linear_lopo_cluster(x[ii],y[ii],ii); 1-sum((y[ii]-pp)^2)/sum((y[ii]-mean(y[ii]))^2)})
    bi<-bi+1L; pred_boot_internal[[bi]]<-data.table(cohort=co,axis=axis,replicate=seq_len(B),Q2=qb)
    for(model in c("LINEAR","HUBER","SPLINE_DF3")) {
      pp <- if(model=="LINEAR")linpred else model_lopo(x,y,model); mm<-pred_metrics(x,y,pp,model)
      ri<-ri+1L; pred_rows[[ri]] <- data.table(cohort=co,axis=axis,model=model,n=length(y),pearson=mm["pearson"],spearman=mm["spearman"],in_sample_R2=mm["in_sample_R2"],LOPO_Q2=mm["Q2"],MAE=mm["MAE"],RMSE=mm["RMSE"],Q2_ci_low=if(model=="LINEAR")quantile(qb,.025,na.rm=TRUE) else NA_real_,Q2_ci_high=if(model=="LINEAR")quantile(qb,.975,na.rm=TRUE) else NA_real_,pearson_permutation_p=if(model=="LINEAR")p_p else NA_real_,spearman_permutation_p=if(model=="LINEAR")p_s else NA_real_,permutation_replicates=if(model=="LINEAR")B else NA_integer_,bootstrap_replicates=if(model=="LINEAR")B else NA_integer_)
    }
  }
  conditioned[[co]] <- rdt
}
pred <- rbindlist(pred_rows); pred[model=="LINEAR",pearson_FDR:=p.adjust(pearson_permutation_p,"BH"),by=cohort]; pred[model=="LINEAR",spearman_FDR:=p.adjust(spearman_permutation_p,"BH"),by=cohort]
write_tab(pred,"02_fingerprint/AXIS_PREDICTABILITY_LOPO.tsv")
write_tab(rbindlist(pred_boot_internal),"02_fingerprint/AXIS_PREDICTABILITY_BOOTSTRAP_INTERNAL.tsv")
cond <- rbindlist(conditioned); setorder(cond,cohort,patient_id); write_tab(cond,"02_fingerprint/MES_CONDITIONED_FINGERPRINT.tsv")
qa_fp <- rbindlist(lapply(c("CGGA325","CGGA693"),function(co) rbindlist(lapply(names(axis_z_cols),function(axis){
  rr<-cond[cohort==co][[paste0("residual_z_",axis)]]; bb<-cond[cohort==co]$bulk_mes_z
  data.table(cohort=co,axis=axis,n=sum(is.finite(rr)),mean=mean(rr,na.rm=TRUE),SD=sd(rr,na.rm=TRUE),pearson_with_MES=cor(rr,bb,use="complete.obs"),spearman_with_MES=cor(rr,bb,use="complete.obs",method="spearman"),min=min(rr,na.rm=TRUE),max=max(rr,na.rm=TRUE),n_abs_gt3=sum(abs(rr)>3,na.rm=TRUE),missing_n=sum(!is.finite(rr)),outlier_action="RETAIN_ALL")
}))))
write_tab(qa_fp,"02_fingerprint/FINGERPRINT_QA.tsv")
multi <- pred[model=="LINEAR",.(mean_LOPO_Q2=mean(LOPO_Q2),median_LOPO_Q2=median(LOPO_Q2),n_Q2_lt_0=sum(LOPO_Q2<0),n_Q2_lt_0_25=sum(LOPO_Q2<.25),n_Q2_lt_0_50=sum(LOPO_Q2<.5),n_axes=.N),by=cohort]
write_tab(multi,"02_fingerprint/MULTIAXIS_PREDICTABILITY_SUMMARY.tsv")

# Pair and nearest-neighbour machinery.
nearest_indices <- function(score,ids,exclude_key=ids) {
  n<-length(score); u<-sort(unique(score)); gi<-match(score,u); groups<-split(seq_len(n),gi); ans<-integer(n)
  for(g in seq_along(groups)) {
    members<-groups[[g]]
    for(i in members) {
      same_score<-members[exclude_key[members]!=exclude_key[i]]
      if(length(same_score)){ans[i]<-same_score[order(ids[same_score])][1];next}
      cand<-integer(); dd<-numeric()
      if(g>1L){cc<-groups[[g-1L]][exclude_key[groups[[g-1L]]]!=exclude_key[i]];cand<-c(cand,cc);dd<-c(dd,rep(score[i]-u[g-1L],length(cc)))}
      if(g<length(groups)){cc<-groups[[g+1L]][exclude_key[groups[[g+1L]]]!=exclude_key[i]];cand<-c(cand,cc);dd<-c(dd,rep(u[g+1L]-score[i],length(cc)))}
      if(!length(cand))stop("No distinct-patient nearest-neighbour candidate")
      md<-min(dd); cc<-cand[dd==md]; ans[i]<-cc[order(ids[cc])][1]
    }
  }
  ans
}
pair_table <- function(ids,score,F,R=NULL,Fenv=NULL,Fraw=NULL,original_ids=ids) {
  ix<-t(combn(seq_along(ids),2));ix<-ix[original_ids[ix[,1]]!=original_ids[ix[,2]],,drop=FALSE]; diffF<-F[ix[,1],,drop=FALSE]-F[ix[,2],,drop=FALSE]
  out<-data.table(patient_i=ids[ix[,1]],patient_j=ids[ix[,2]],delta_score=abs(score[ix[,1]]-score[ix[,2]]),distance_euclidean_raw=sqrt(rowSums(diffF^2)),distance_manhattan_raw=rowSums(abs(diffF)))
  if(!identical(original_ids,ids))out[,`:=`(original_patient_i=original_ids[ix[,1]],original_patient_j=original_ids[ix[,2]])]
  if(!is.null(R)){dr<-R[ix[,1],,drop=FALSE]-R[ix[,2],,drop=FALSE];out[,distance_euclidean_residual:=sqrt(rowSums(dr^2))]}
  if(!is.null(Fenv)){de<-Fenv[ix[,1],,drop=FALSE]-Fenv[ix[,2],,drop=FALSE];out[,distance_environment_only:=sqrt(rowSums(de^2))]}
  if(!is.null(Fraw)){du<-Fraw[ix[,1],,drop=FALSE]-Fraw[ix[,2],,drop=FALSE];out[,distance_raw_ecology_unpruned:=sqrt(rowSums(du^2))]}
  out
}
nn_table <- function(ids,score,F,R=NULL,Fenv=NULL,Fraw=NULL,original_ids=ids) {
  j<-nearest_indices(score,ids,original_ids); df<-F-F[j,,drop=FALSE]
  out<-data.table(patient_id=ids,nearest_patient_id=ids[j],delta_score=abs(score-score[j]),distance_euclidean_raw=sqrt(rowSums(df^2)),distance_manhattan_raw=rowSums(abs(df)))
  if(!is.null(R)){dr<-R-R[j,,drop=FALSE];out[,distance_euclidean_residual:=sqrt(rowSums(dr^2))]}
  if(!is.null(Fenv)){de<-Fenv-Fenv[j,,drop=FALSE];out[,distance_environment_only:=sqrt(rowSums(de^2))]}
  if(!is.null(Fraw)){du<-Fraw-Fraw[j,,drop=FALSE];out[,distance_raw_ecology_unpruned:=sqrt(rowSums(du^2))]}
  out
}
make_primary_mats <- function(x,key=x$patient_id) {
  F<-as.matrix(x[,.(malignant=zv(malignant_mes_raw),myeloid=zv(myeloid_fraction),hypoxia=zv(hypoxia_pruned),matrix=zv(matrix_pruned),vascular=zv(vascular_stromal_fraction))])
  Fenv<-F[,c("myeloid","hypoxia","matrix","vascular"),drop=FALSE]
  Fraw<-as.matrix(x[,.(malignant=zv(malignant_mes_raw),myeloid=zv(myeloid_fraction),hypoxia=zv(hypoxia_raw),matrix=zv(matrix_raw),vascular=zv(vascular_stromal_fraction))])
  xscore<-zv(x$bulk_mes_raw); R<-sapply(seq_len(ncol(F)),function(k)zv(F[,k]-linear_lopo_cluster(xscore,F[,k],key)));colnames(R)<-colnames(F)
  list(score=xscore,F=F,Fenv=Fenv,Fraw=Fraw,R=R)
}
all_pairs<-list(); nn_all<-list(); threshold_rows<-list(); hero_rows<-list(); perm_rows<-list(); perm_summary<-list(); boot_rep<-list(); tri<-0L; pri<-0L;br<-0L
for(coi in seq_along(c("CGGA325","CGGA693"))) {
  co<-c("CGGA325","CGGA693")[coi]; x<-copy(d[cohort==co]); m<-make_primary_mats(x)
  pp<-pair_table(x$patient_id,m$score,m$F,m$R,m$Fenv,m$Fraw);pp[,cohort:=co];setcolorder(pp,c("cohort",setdiff(names(pp),"cohort")));all_pairs[[co]]<-pp
  nn<-nn_table(x$patient_id,m$score,m$F,m$R,m$Fenv,m$Fraw);nn[,`:=`(cohort=co,bulk_mes_z=m$score,nearest_bulk_mes_z=m$score[match(nearest_patient_id,x$patient_id)])];setcolorder(nn,c("cohort","patient_id","nearest_patient_id","bulk_mes_z","nearest_bulk_mes_z",setdiff(names(nn),c("cohort","patient_id","nearest_patient_id","bulk_mes_z","nearest_bulk_mes_z"))));nn_all[[co]]<-nn
  allmed<-median(pp$distance_euclidean_raw)
  for(th in THRESHOLDS) {
    z<-pp[delta_score<=th]; tri<-tri+1L; q<-qtiles(z$distance_euclidean_raw)
    threshold_rows[[tri]]<-data.table(cohort=co,threshold_sd=th,is_primary=th==.10,n_pairs=nrow(z),n_unique_patients=uniqueN(c(z$patient_i,z$patient_j)),percent_cohort_represented=100*uniqueN(c(z$patient_i,z$patient_j))/nrow(x),median_delta_score=median(z$delta_score),median_configuration_distance=q[2],Q1_configuration_distance=q[1],Q3_configuration_distance=q[3],P90_configuration_distance=q[4],max_configuration_distance=max(z$distance_euclidean_raw),median_environment_only_distance=median(z$distance_environment_only),median_residual_fingerprint_distance=median(z$distance_euclidean_residual),median_manhattan_distance=median(z$distance_manhattan_raw),median_raw_ecology_unpruned_distance=median(z$distance_raw_ecology_unpruned),all_pair_median_configuration_distance=allmed,divergent_fraction_gt_all_pair_median=mean(z$distance_euclidean_raw>allmed))
  }
  cand<-pp[delta_score<=.10][order(-distance_euclidean_raw,delta_score,patient_i,patient_j)][1]
  a<-x[patient_id==cand$patient_i]; b2<-x[patient_id==cand$patient_j]
  hero_rows[[co]]<-data.table(cohort=co,patient_i=cand$patient_i,patient_j=cand$patient_j,bulk_mes_raw_i=a$bulk_mes_raw,bulk_mes_raw_j=b2$bulk_mes_raw,bulk_mes_z_i=a$bulk_mes_z,bulk_mes_z_j=b2$bulk_mes_z,delta_mes=cand$delta_score,
    malignant_mes_raw_i=a$malignant_mes_raw,malignant_mes_raw_j=b2$malignant_mes_raw,malignant_mes_z_i=a$malignant_mes_z,malignant_mes_z_j=b2$malignant_mes_z,
    myeloid_fraction_i=a$myeloid_fraction,myeloid_fraction_j=b2$myeloid_fraction,myeloid_z_i=a$myeloid_z,myeloid_z_j=b2$myeloid_z,
    hypoxia_pruned_i=a$hypoxia_pruned,hypoxia_pruned_j=b2$hypoxia_pruned,hypoxia_z_i=a$hypoxia_z,hypoxia_z_j=b2$hypoxia_z,
    matrix_pruned_i=a$matrix_pruned,matrix_pruned_j=b2$matrix_pruned,matrix_z_i=a$matrix_z,matrix_z_j=b2$matrix_z,
    vascular_stromal_fraction_i=a$vascular_stromal_fraction,vascular_stromal_fraction_j=b2$vascular_stromal_fraction,vascular_z_i=a$vascular_z,vascular_z_j=b2$vascular_z,
    raw_fingerprint_distance=cand$distance_euclidean_raw,residual_fingerprint_distance=cand$distance_euclidean_residual,environment_only_distance=cand$distance_environment_only,selection_rule="MAX_RAW_DISTANCE_WITHIN_DELTA_MES_LE_0.10; TIE_SMALLER_DELTA_THEN_LEXICOGRAPHIC")
  set.seed(SEEDS$nn_perm+coi)
  pmed<-numeric(B);penv<-numeric(B)
  for(k in seq_len(B)){sp<-sample(m$score);j<-nearest_indices(sp,x$patient_id);pmed[k]<-median(sqrt(rowSums((m$F-m$F[j,,drop=FALSE])^2)));penv[k]<-median(sqrt(rowSums((m$Fenv-m$Fenv[j,,drop=FALSE])^2)))}
  perm_rows[[co]]<-data.table(cohort=co,permutation=seq_len(B),median_nn_configuration_distance=pmed,median_nn_environment_distance=penv)
  obs<-median(nn$distance_euclidean_raw);obse<-median(nn$distance_environment_only)
  perm_summary[[co]]<-data.table(cohort=co,n_patients=nrow(x),observed_median_nn_distance=obs,null_median_of_medians=median(pmed),retention_ratio=obs/median(pmed),one_minus_retention=1-obs/median(pmed),finite_permutation_p=(1+sum(pmed<=obs))/(B+1),observed_median_nn_environment=obse,null_median_environment=median(penv),environment_retention_ratio=obse/median(penv),environment_permutation_p=(1+sum(penv<=obse))/(B+1),permutations=B)
  set.seed(SEEDS$iso_boot+coi)
  for(k in seq_len(B)) {
    ii<-sample.int(nrow(x),nrow(x),TRUE); xb<-x[ii]; ids<-sprintf("%s__BOOT%03d",xb$patient_id,seq_len(nrow(xb))); mb<-make_primary_mats(xb,xb$patient_id)
    pb<-pair_table(ids,mb$score,mb$F,mb$R,mb$Fenv,mb$Fraw,xb$patient_id);nb<-nn_table(ids,mb$score,mb$F,mb$R,mb$Fenv,mb$Fraw,xb$patient_id);amed<-median(pb$distance_euclidean_raw)
    br<-br+1L;boot_rep[[br]]<-data.table(cohort=co,replicate=k,threshold_sd=NA_real_,metric=c("median_nn_distance","median_nn_environment","median_nearest_delta_score"),value=c(median(nb$distance_euclidean_raw),median(nb$distance_environment_only),median(nb$delta_score)))
    for(th in THRESHOLDS){zz<-pb[delta_score<=th];br<-br+1L;boot_rep[[br]]<-data.table(cohort=co,replicate=k,threshold_sd=th,metric=c("iso_n_pairs","iso_percent_patients","iso_median_distance","iso_P90_distance","iso_divergent_fraction","iso_median_environment","iso_median_residual","iso_median_manhattan","iso_median_raw_ecology_unpruned"),value=c(nrow(zz),100*uniqueN(c(zz$original_patient_i,zz$original_patient_j))/uniqueN(xb$patient_id),median(zz$distance_euclidean_raw),quantile(zz$distance_euclidean_raw,.9,type=7),mean(zz$distance_euclidean_raw>amed),median(zz$distance_environment_only),median(zz$distance_euclidean_residual),median(zz$distance_manhattan_raw),median(zz$distance_raw_ecology_unpruned)))}
  }
}
pairs<-rbindlist(all_pairs);setnames(pairs,"delta_score","delta_mes");write_tab(pairs,"03_iso_score/ALL_PATIENT_PAIR_DISTANCES.tsv")
nnmes<-rbindlist(nn_all);setnames(nnmes,"delta_score","delta_mes");write_tab(nnmes,"03_iso_score/NEAREST_MES_NEIGHBOURS.tsv")
thsum<-rbindlist(threshold_rows);write_tab(thsum,"03_iso_score/ISO_MES_THRESHOLD_SUMMARY.tsv")
permd<-rbindlist(perm_rows);write_tab(permd,"03_iso_score/NEAREST_MES_PERMUTATION.tsv")
perms<-rbindlist(perm_summary);write_tab(perms,"03_iso_score/NEAREST_MES_PERMUTATION_SUMMARY.tsv")
heros<-rbindlist(hero_rows);write_tab(heros,"03_iso_score/ISO_MES_HERO_PAIRS.tsv")
brep<-rbindlist(boot_rep);write_tab(brep,"03_iso_score/PATIENT_BOOTSTRAP_REPLICATES_INTERNAL.tsv")
bsummary<-brep[,.(effect=median(value,na.rm=TRUE),ci_low=quantile(value,.025,na.rm=TRUE,type=7),ci_high=quantile(value,.975,na.rm=TRUE,type=7),bootstrap_replicates=.N),by=.(cohort,threshold_sd,metric)]
# The point estimate is always the observed-cohort statistic; bootstrap medians
# are not substituted for the observed effect.
metric_map<-c(iso_n_pairs="n_pairs",iso_percent_patients="percent_cohort_represented",iso_median_distance="median_configuration_distance",iso_P90_distance="P90_configuration_distance",iso_divergent_fraction="divergent_fraction_gt_all_pair_median",iso_median_environment="median_environment_only_distance",iso_median_residual="median_residual_fingerprint_distance",iso_median_manhattan="median_manhattan_distance",iso_median_raw_ecology_unpruned="median_raw_ecology_unpruned_distance")
obs_boot_points<-rbindlist(list(
  nnmes[,.(effect=median(distance_euclidean_raw)),by=cohort][,`:=`(threshold_sd=NA_real_,metric="median_nn_distance")],
  nnmes[,.(effect=median(distance_environment_only)),by=cohort][,`:=`(threshold_sd=NA_real_,metric="median_nn_environment")],
  nnmes[,.(effect=median(delta_mes)),by=cohort][,`:=`(threshold_sd=NA_real_,metric="median_nearest_delta_score")],
  rbindlist(lapply(names(metric_map),function(mm)thsum[,.(cohort,threshold_sd,metric=mm,effect=get(metric_map[[mm]]))]))
),fill=TRUE)
bsummary[obs_boot_points,on=.(cohort,threshold_sd,metric),effect:=i.effect]
# Add conditional retention-ratio intervals using the frozen 5,000-permutation denominator.
for(co in c("CGGA325","CGGA693")){den<-perms[cohort==co,null_median_of_medians];denv<-perms[cohort==co,null_median_environment];z<-brep[cohort==co&metric=="median_nn_distance",value/den];ze<-brep[cohort==co&metric=="median_nn_environment",value/denv];bsummary<-rbind(bsummary,data.table(cohort=co,threshold_sd=NA_real_,metric=c("retention_ratio_conditional_on_permutation_null","environment_retention_ratio_conditional_on_permutation_null"),effect=c(perms[cohort==co,retention_ratio],perms[cohort==co,environment_retention_ratio]),ci_low=c(quantile(z,.025),quantile(ze,.025)),ci_high=c(quantile(z,.975),quantile(ze,.975)),bootstrap_replicates=B))}
write_tab(bsummary,"03_iso_score/PATIENT_BOOTSTRAP_SUMMARY.tsv")

# Exact v1.4 additive provenance link and hero additions.
prov<-fread(file.path(V14,"02_additive_contributions/STRICT_BULK_CLASS_CONTRIBUTIONS.tsv"))[cohort%in%c("CGGA325","CGGA693")]
prov<-prov[,.(cohort,patient_id,C_MALIGNANT,C_ECOLOGICAL,C_SHARED,C_UNSTABLE,MES_ADDITIVE)]
prov<-merge(prov,d[,.(cohort,patient_id,bulk_mes_raw,bulk_mes_z,mes_q4)],by=c("cohort","patient_id"),all.y=TRUE)
prov[,additivity_error:=abs(C_MALIGNANT+C_ECOLOGICAL+C_SHARED+C_UNSTABLE-bulk_mes_raw)]
if(nrow(prov)!=181L||max(prov$additivity_error)>1e-10)stop("v1.4 provenance additivity failure")
pmh<-prov[mes_q4==TRUE];setorder(pmh,cohort,bulk_mes_z,patient_id);pmh[,mes_high_order:=seq_len(.N),by=cohort]
write_tab(pmh,"05_provenance_link/MES_HIGH_PROVENANCE_REALIZATION.tsv")
for(co in c("CGGA325","CGGA693")){
  for(side in c("i","j")){
    ids<-heros[cohort==co][[paste0("patient_",side)]];pr<-prov[cohort==co&patient_id==ids]
    for(nm in c("C_MALIGNANT","C_ECOLOGICAL","C_SHARED","C_UNSTABLE"))heros[cohort==co,(paste0(nm,"_",side)):=pr[[nm]]]
  }
}
write_tab(heros,"03_iso_score/ISO_MES_HERO_PAIRS.tsv")

# State-specific fingerprints and controls.
state_def <- list(
  CANONICAL_MES=list(bulk="bulk_mes_raw",mal="malignant_mes_raw",cov="coverage_fraction_CANONICAL_MES"),
  AC_LIKE=list(bulk="bulk_AC_LIKE",mal="malignant_AC_LIKE",cov="coverage_fraction_AC_LIKE"),
  OPC_LIKE=list(bulk="bulk_OPC_LIKE",mal="malignant_OPC_LIKE",cov="coverage_fraction_OPC_LIKE"),
  NPC_LIKE=list(bulk="bulk_NPC_LIKE",mal="malignant_NPC_LIKE",cov="coverage_fraction_NPC_LIKE"),
  PROLIFERATION=list(bulk="bulk_PROLIFERATION",mal="malignant_PROLIFERATION",cov="coverage_fraction_PROLIFERATION")
)
state_mats <- function(x,st) {
  def<-state_def[[st]]; score<-zv(x[[def$bulk]])
  F<-as.matrix(data.table(malignant=zv(x[[def$mal]]),myeloid=zv(x$myeloid_fraction),hypoxia=zv(x[[paste0("hypoxia_pruned_",st)]]),matrix=zv(x[[paste0("matrix_pruned_",st)]]),vascular=zv(x$vascular_stromal_fraction)))
  Fenv<-F[,c("myeloid","hypoxia","matrix","vascular"),drop=FALSE]
  Fraw<-as.matrix(data.table(malignant=zv(x[[def$mal]]),myeloid=zv(x$myeloid_fraction),hypoxia=zv(x$hypoxia_raw),matrix=zv(x$matrix_raw),vascular=zv(x$vascular_stromal_fraction)))
  list(score=score,F=F,Fenv=Fenv,Fraw=Fraw)
}
state_nn_summary<-list();state_nn_pat<-list();state_perm_rep<-list();state_iso<-list();state_pred<-list();sri<-0L;spri<-0L;sii<-0L;sapi<-0L
state_null <- list()
for(coi in seq_along(c("CGGA325","CGGA693"))) {
  co<-c("CGGA325","CGGA693")[coi];x<-copy(d[cohort==co])
  for(si in seq_along(state_def)) {
    st<-names(state_def)[si];def<-state_def[[st]]; evaluable<-all(is.finite(x[[def$mal]]))&&min(x[[def$cov]],na.rm=TRUE)>=.5
    status<-if(evaluable)"EVALUATED" else "NOT_EVALUABLE_LOW_GENE_COVERAGE"
    if(!evaluable){sri<-sri+1L;state_nn_summary[[sri]]<-data.table(cohort=co,state=st,status=status,n=nrow(x),median_nn_distance=NA_real_,Q1=NA_real_,Q3=NA_real_,P90=NA_real_,max=NA_real_,median_delta_score=NA_real_,null_median=NA_real_,retention_ratio=NA_real_,one_minus_retention=NA_real_,permutation_p=NA_real_,permutations=0L);next}
    m<-state_mats(x,st);nn<-nn_table(x$patient_id,m$score,m$F,NULL,m$Fenv,m$Fraw);nn[,`:=`(cohort=co,state=st,status=status,score_z=m$score,nearest_score_z=m$score[match(nearest_patient_id,x$patient_id)])];spri<-spri+1L;state_nn_pat[[spri]]<-nn
    set.seed(SEEDS$state_perm+coi*100L+si);pmed<-numeric(B);penv<-numeric(B)
    for(k in seq_len(B)){sp<-sample(m$score);j<-nearest_indices(sp,x$patient_id);pmed[k]<-median(sqrt(rowSums((m$F-m$F[j,,drop=FALSE])^2)));penv[k]<-median(sqrt(rowSums((m$Fenv-m$Fenv[j,,drop=FALSE])^2)))}
    sapi<-sapi+1L;state_perm_rep[[sapi]]<-data.table(cohort=co,state=st,replicate=seq_len(B),median_nn_distance=pmed,median_nn_environment=penv)
    state_null[[paste(co,st,sep="::")]]<-c(raw=median(pmed),env=median(penv))
    q<-qtiles(nn$distance_euclidean_raw);sri<-sri+1L;state_nn_summary[[sri]]<-data.table(cohort=co,state=st,status=status,n=nrow(x),median_nn_distance=q[2],Q1=q[1],Q3=q[3],P90=q[4],max=max(nn$distance_euclidean_raw),median_delta_score=median(nn$delta_score),null_median=median(pmed),retention_ratio=q[2]/median(pmed),one_minus_retention=1-q[2]/median(pmed),permutation_p=(1+sum(pmed<=q[2]))/(B+1),permutations=B,median_nn_environment=median(nn$distance_environment_only),environment_null_median=median(penv),environment_retention_ratio=median(nn$distance_environment_only)/median(penv),environment_permutation_p=(1+sum(penv<=median(nn$distance_environment_only)))/(B+1))
    pp<-pair_table(x$patient_id,m$score,m$F,NULL,m$Fenv,m$Fraw);amed<-median(pp$distance_euclidean_raw)
    for(th in THRESHOLDS){zz<-pp[delta_score<=th];sii<-sii+1L;state_iso[[sii]]<-data.table(cohort=co,state=st,status=status,threshold_sd=th,n_pairs=nrow(zz),n_unique_patients=uniqueN(c(zz$patient_i,zz$patient_j)),median_distance=median(zz$distance_euclidean_raw),P90_distance=quantile(zz$distance_euclidean_raw,.9,type=7),divergent_fraction=mean(zz$distance_euclidean_raw>amed),median_environment=median(zz$distance_environment_only),median_raw_ecology_unpruned=median(zz$distance_raw_ecology_unpruned))}
    for(ai in seq_len(ncol(m$F))){axis<-colnames(m$F)[ai];y<-m$F[,ai];set.seed(SEEDS$axis_boot+coi*10000L+si*100L+ai);qb<-replicate(B,{ii<-sample.int(nrow(x),nrow(x),TRUE);ppp<-linear_lopo_cluster(m$score[ii],y[ii],ii);1-sum((y[ii]-ppp)^2)/sum((y[ii]-mean(y[ii]))^2)});for(model in c("LINEAR","HUBER","SPLINE_DF3")){pr<-model_lopo(m$score,y,model);mm<-pred_metrics(m$score,y,pr,model);state_pred[[length(state_pred)+1L]]<-data.table(cohort=co,state=st,axis=axis,model=model,n=nrow(x),LOPO_Q2=mm["Q2"],MAE=mm["MAE"],RMSE=mm["RMSE"],pearson=mm["pearson"],spearman=mm["spearman"],in_sample_R2=mm["in_sample_R2"],Q2_ci_low=if(model=="LINEAR")quantile(qb,.025,na.rm=TRUE)else NA_real_,Q2_ci_high=if(model=="LINEAR")quantile(qb,.975,na.rm=TRUE)else NA_real_,bootstrap_replicates=if(model=="LINEAR")B else NA_integer_)}}
  }
}
state_sum<-rbindlist(state_nn_summary,fill=TRUE);write_tab(state_sum,"04_state_specificity/STATE_SPECIFICITY_NEAREST_NEIGHBOUR.tsv")
write_tab(rbindlist(state_nn_pat,fill=TRUE),"04_state_specificity/STATE_SPECIFICITY_NEAREST_NEIGHBOUR_PATIENTS_INTERNAL.tsv")
write_tab(rbindlist(state_perm_rep),"04_state_specificity/STATE_SPECIFICITY_PERMUTATION_INTERNAL.tsv")
write_tab(rbindlist(state_iso),"04_state_specificity/STATE_SPECIFICITY_ISO_THRESHOLD_SUMMARY.tsv")
write_tab(rbindlist(state_pred),"04_state_specificity/STATE_SPECIFICITY_AXIS_PREDICTABILITY.tsv")

# Paired patient bootstrap for state specificity.
state_boot_rows<-list();sbr<-0L
for(coi in seq_along(c("CGGA325","CGGA693"))) {
  co<-c("CGGA325","CGGA693")[coi];x<-copy(d[cohort==co]);set.seed(SEEDS$state_boot+coi)
  eval_states<-names(state_def)[vapply(names(state_def),function(st){def<-state_def[[st]];all(is.finite(x[[def$mal]]))&&min(x[[def$cov]],na.rm=TRUE)>=.5},logical(1))]
  for(k in seq_len(B)) {
    ii<-sample.int(nrow(x),nrow(x),TRUE);xb<-x[ii];ids<-sprintf("%s__BOOT%03d",xb$patient_id,seq_len(nrow(xb)));vals<-list()
    for(st in eval_states){m<-state_mats(xb,st);nn<-nn_table(ids,m$score,m$F,NULL,m$Fenv,m$Fraw,xb$patient_id);vals[[st]]<-c(raw=median(nn$distance_euclidean_raw),env=median(nn$distance_environment_only))}
    mes<-vals[["CANONICAL_MES"]]
    for(st in eval_states){sbr<-sbr+1L;state_boot_rows[[sbr]]<-data.table(cohort=co,state=st,replicate=k,median_nn_distance=vals[[st]]["raw"],median_nn_environment=vals[[st]]["env"],delta_mes_minus_control=if(st=="CANONICAL_MES")0 else mes["raw"]-vals[[st]]["raw"],delta_mes_minus_control_environment=if(st=="CANONICAL_MES")0 else mes["env"]-vals[[st]]["env"],retention_ratio_conditional=vals[[st]]["raw"]/state_null[[paste(co,st,sep="::")]]["raw"],environment_retention_ratio_conditional=vals[[st]]["env"]/state_null[[paste(co,st,sep="::")]]["env"])}
  }
}
sb<-rbindlist(state_boot_rows)
sb[,mes_retention_boot:=retention_ratio_conditional[state=="CANONICAL_MES"],by=.(cohort,replicate)]
sb[,mes_environment_retention_boot:=environment_retention_ratio_conditional[state=="CANONICAL_MES"],by=.(cohort,replicate)]
sb[,`:=`(delta_retention_mes_minus_control=mes_retention_boot-retention_ratio_conditional,delta_environment_retention_mes_minus_control=mes_environment_retention_boot-environment_retention_ratio_conditional)]
write_tab(sb,"04_state_specificity/STATE_SPECIFICITY_BOOTSTRAP_REPLICATES_INTERNAL.tsv")
boot_summary <- function(value_col) sb[,.(bootstrap_median=median(get(value_col),na.rm=TRUE),ci_low=quantile(get(value_col),.025,na.rm=TRUE),ci_high=quantile(get(value_col),.975,na.rm=TRUE),bootstrap_replicates=.N),by=.(cohort,state)]
point_raw<-state_sum[,.(cohort,state,status,effect=median_nn_distance)];point_env<-state_sum[,.(cohort,state,status,effect=median_nn_environment)]
point_rr<-state_sum[,.(cohort,state,status,effect=retention_ratio)];point_err<-state_sum[,.(cohort,state,status,effect=environment_retention_ratio)]
mesraw<-state_sum[state=="CANONICAL_MES",.(cohort,mes_raw=median_nn_distance,mes_env=median_nn_environment,mes_rr=retention_ratio,mes_err=environment_retention_ratio)]
delta_point<-merge(state_sum[,.(cohort,state,status,control_raw=median_nn_distance,control_env=median_nn_environment,control_rr=retention_ratio,control_err=environment_retention_ratio)],mesraw,by="cohort",all.x=TRUE)
raw1<-merge(point_raw,boot_summary("median_nn_distance"),by=c("cohort","state"),all.x=TRUE);raw1[,metric:="median_nn_distance"]
raw2<-merge(delta_point[,.(cohort,state,status,effect=mes_raw-control_raw)],boot_summary("delta_mes_minus_control"),by=c("cohort","state"),all.x=TRUE);raw2[,metric:="delta_mes_minus_control"]
raw3<-merge(point_rr,boot_summary("retention_ratio_conditional"),by=c("cohort","state"),all.x=TRUE);raw3[,metric:="retention_ratio"]
raw4<-merge(delta_point[,.(cohort,state,status,effect=mes_rr-control_rr)],boot_summary("delta_retention_mes_minus_control"),by=c("cohort","state"),all.x=TRUE);raw4[,metric:="delta_retention_mes_minus_control"]
sout<-rbindlist(list(raw1,raw2,raw3,raw4),fill=TRUE,use.names=TRUE);setcolorder(sout,c("metric","cohort","state","status","effect","bootstrap_median","ci_low","ci_high","bootstrap_replicates"));write_tab(sout,"04_state_specificity/STATE_SPECIFICITY_BOOTSTRAP.tsv")
env1<-merge(point_env,boot_summary("median_nn_environment"),by=c("cohort","state"),all.x=TRUE);env1[,metric:="median_nn_environment"]
env2<-merge(delta_point[,.(cohort,state,status,effect=mes_env-control_env)],boot_summary("delta_mes_minus_control_environment"),by=c("cohort","state"),all.x=TRUE);env2[,metric:="delta_mes_minus_control_environment"]
env3<-merge(point_err,boot_summary("environment_retention_ratio_conditional"),by=c("cohort","state"),all.x=TRUE);env3[,metric:="environment_retention_ratio"]
env4<-merge(delta_point[,.(cohort,state,status,effect=mes_err-control_err)],boot_summary("delta_environment_retention_mes_minus_control"),by=c("cohort","state"),all.x=TRUE);env4[,metric:="delta_environment_retention_mes_minus_control"]
envout<-rbindlist(list(env1,env2,env3,env4),fill=TRUE,use.names=TRUE);setcolorder(envout,c("metric","cohort","state","status","effect","bootstrap_median","ci_low","ci_high","bootstrap_replicates"));write_tab(envout,"04_state_specificity/STATE_SPECIFICITY_ENVIRONMENT_ONLY.tsv")

# Existing pseudo-bulk controlled support.
pbdesign<-fread(file.path(V11,"04_pseudobulk_benchmark/PSEUDOBULK_DESIGN.tsv"));pbresp<-fread(file.path(V11,"04_pseudobulk_benchmark/PSEUDOBULK_SCORE_RESPONSES.tsv"))[score_name=="canonical_MES"]
pb<-merge(pbdesign,pbresp[,.(pseudobulk_id,canonical_MES=score)],by="pseudobulk_id",all.x=TRUE)
metadata_ok<-nrow(pb)==8100L&&all(c("patient_key","composition_fraction","compartment","canonical_MES")%in%names(pb))&&!anyNA(pb[,.(patient_key,composition_fraction,compartment,canonical_MES)])
pbpairs<-list();pbi<-0L
if(metadata_ok){
  pa<-pb[experiment=="A_FIXED_MALIGNANT"]
  pa[,`:=`(fixed_state_group=paste(study,patient_key,sep="::"),malignant_fraction=1-composition_fraction,myeloid_fraction=ifelse(compartment=="myeloid",composition_fraction,0),vascular_stromal_fraction=ifelse(compartment%in%c("endothelial","pericyte"),composition_fraction,0),other_fraction=ifelse(compartment%in%c("lymphoid","oligodendrocyte_other"),composition_fraction,0))]
  for(gr in unique(pa$fixed_state_group)){z<-pa[fixed_state_group==gr];ix<-t(combn(seq_len(nrow(z)),2));delta<-abs(z$canonical_MES[ix[,1]]-z$canonical_MES[ix[,2]]);cm<-as.matrix(z[,.(malignant_fraction,myeloid_fraction,vascular_stromal_fraction,other_fraction)]);cd<-sqrt(rowSums((cm[ix[,1],,drop=FALSE]-cm[ix[,2],,drop=FALSE])^2));keep<-delta<=.10&cd>0;if(any(keep)){pbi<-pbi+1L;pbpairs[[pbi]]<-data.table(fixed_state_group=gr,mixture_i=z$pseudobulk_id[ix[keep,1]],mixture_j=z$pseudobulk_id[ix[keep,2]],compartment=z$compartment[ix[keep,1]],delta_mes=delta[keep],composition_distance=cd[keep],composition_fraction_i=z$composition_fraction[ix[keep,1]],composition_fraction_j=z$composition_fraction[ix[keep,2]],canonical_mes_i=z$canonical_MES[ix[keep,1]],canonical_mes_j=z$canonical_MES[ix[keep,2]],status="EVALUATED_EXISTING_REAL_COUNT_MIXTURES",interpretation_boundary="CONTROLLED_SUPPORT_NOT_SIMULATED_PATIENT_EVIDENCE")}}
}
if(metadata_ok&&length(pbpairs)) write_tab(rbindlist(pbpairs),"07_pseudobulk_control/PSEUDOBULK_ISO_MES_CONFIGURATION.tsv") else write_tab(data.table(status="NOT_EVALUABLE_EXISTING_MIXTURE_METADATA",reason="Required identity/composition/response metadata incomplete"),"07_pseudobulk_control/PSEUDOBULK_ISO_MES_CONFIGURATION.tsv")
write_tab(data.table(check=c("rows_8100","fixed_state_identity","composition_metadata","canonical_mes_response","mixtures_regenerated"),status=c(ifelse(nrow(pb)==8100,"PASS","FAIL"),ifelse(all(nzchar(pb$patient_key)),"PASS","FAIL"),ifelse(all(is.finite(pb$composition_fraction)),"PASS","FAIL"),ifelse(all(is.finite(pb$canonical_MES)),"PASS","FAIL"),"PASS_NO_REGENERATION"),detail=c(nrow(pb),uniqueN(pb$patient_key),"existing design fractions","existing v1.1 response","read-only reuse")),"07_pseudobulk_control/PSEUDOBULK_CONTROL_QA.tsv")

writeLines(capture.output(sessionInfo()),file.path(RUN,"session_info/CORE_ANALYSIS_SESSION_INFO.txt"))
cat("FORMAL_ANALYSIS_COMPLETE",format(Sys.time(),tz="UTC",usetz=TRUE),"\n")
