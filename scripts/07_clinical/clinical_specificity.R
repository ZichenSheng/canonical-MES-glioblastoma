#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
suppressPackageStartupMessages({
  library(BayesPrism); library(Matrix); library(data.table); library(readxl)
  library(survival); library(MASS)
})

DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset = "data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset = "results")
REPO_ROOT <- Sys.getenv("GBM_MES_REPO_ROOT", unset = ".")
RUN <- file.path(RESULT_ROOT, "clinical_specificity")
OUT <- file.path(RUN, "outputs")
FORMAL <- file.path(RESULT_ROOT, "clinical_icb")
ICB <- file.path(DATA_ROOT, "clinical", "icb")
FIT <- file.path(DATA_ROOT, "prepared", "clinical", "icb", "BAYESPRISM_FULL_ICB_40_FIT.rds")
SIGREG <- file.path(REPO_ROOT, "resources", "signatures", "signature_gene_sets.tsv")
GMT <- file.path(DATA_ROOT, "resources", "msigdb", "h.all.v2026.1.Hs.symbols.gmt")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(RUN, "logs"), recursive=TRUE, showWarnings=FALSE)
w <- function(x,n) fwrite(x,file.path(OUT,n),sep="\t",quote=FALSE,na="NA")
znum <- function(x) {z <- as.numeric(scale(x)); z[!is.finite(z)] <- 0; z}
safe_coef <- function(f, term) if (!is.null(f) && term %in% names(coef(f)) && is.finite(coef(f)[term])) unname(coef(f)[term]) else NA_real_
safe_ph <- function(f, term=NULL) {
  z <- tryCatch(suppressWarnings(cox.zph(f)), error=function(e) NULL)
  if (is.null(z)) return(NA_real_)
  if (is.null(term)) return(unname(z$table["GLOBAL","p"]))
  if (!term %in% rownames(z$table)) return(NA_real_)
  unname(z$table[term,"p"])
}
calc_c <- function(dat, lp) tryCatch(unname(concordance(Surv(time_years,event)~lp, data=dat, reverse=TRUE)$concordance), error=function(e) NA_real_)

# Formal cohort and frozen decomposition.
d0 <- fread(file.path(FORMAL,"02_full_icb_projection/FULL_ICB_PATIENT_COMPONENT_SCORES.tsv"))
stopifnot(nrow(d0)==40L, uniqueN(d0$participant_id)==40L, sum(d0$event)==38L)
raw <- as.data.frame(read_excel(file.path(ICB,"rnaseq_counts_tpms.xlsx"),sheet="sTable6a_rnaseq_counts",skip=1),check.names=FALSE)
names(raw)[1:2] <- c("ensembl","symbol")
sample_cols <- setdiff(names(raw),c("ensembl","symbol"))
lib <- colSums(as.matrix(raw[,sample_cols,drop=FALSE])); names(lib) <- sample_cols
raw$symbol <- toupper(trimws(raw$symbol)); raw <- raw[!is.na(raw$symbol)&nzchar(raw$symbol),,drop=FALSE]
cnt <- as.matrix(raw[,sample_cols,drop=FALSE]); storage.mode(cnt) <- "numeric"; rownames(cnt) <- raw$symbol
cnt <- rowsum(cnt,rownames(cnt),reorder=FALSE,na.rm=TRUE); cnt <- t(cnt)
bulk <- sweep(cnt[d0$sample_id,,drop=FALSE],1,lib[d0$sample_id]/1e6,"/")
fit <- readRDS(FIT)
mal <- get.exp(fit,"type","malignant")[d0$sample_id,,drop=FALSE]
Z <- fit@posterior.initial.cellType@Z
nonmal <- setdiff(dimnames(Z)[[3]],"malignant")
eco <- apply(Z[,,nonmal,drop=FALSE],c(1,2),sum)[d0$sample_id,,drop=FALSE]
theta <- if(!is.null(fit@posterior.theta_f)) fit@posterior.theta_f@theta else fit@posterior.initial.cellType@theta
stopifnot(identical(rownames(bulk),rownames(mal)), identical(rownames(mal),rownames(eco)))

# Certified classifier definitions. The Couturier top-50 rule is inherited verbatim
# from the official project notebook; no current-cohort outcome enters gene selection.
reg <- fread(SIGREG); reg[,gene:=toupper(trimws(gene_symbol_clean))]
canonical <- unique(reg[signature_id=="NEFTEL_MES_LIKE" & !is.na(gene)&nzchar(gene),gene])
ver <- readRDS(file.path(ICB,"GBM_Verhaak_signatures_df.rds"))
neftel <- as.data.table(read_excel(file.path(ICB,"Neftel_et_al_signatures.xlsx"),skip=4))
garofano <- as.data.table(read_excel(file.path(ICB,"Garofano_et_al_signatures.xlsx"),sheet="Tab 26 - Supplementary Table 6j",skip=2))
nomura <- as.data.table(read_excel(file.path(ICB,"Nomura_et_al_signatures.xlsx"),sheet="TableS2",skip=4))
cout <- fread(file.path(ICB,"Couturier_et_al_signatures.csv"))
clean_genes <- function(x) unique(toupper(trimws(as.character(x[!is.na(x)&nzchar(as.character(x))]))))
sets <- list(
  CANONICAL_NEFTEL_MES_LIKE_95_ANCHOR=clean_genes(canonical),
  VERHAAK_MES=clean_genes(ver$gene_name[ver$sig_name=="Mesenchymal"]),
  NEFTEL_MES1=clean_genes(neftel$MES1),
  NEFTEL_MES2=clean_genes(neftel$MES2),
  COUTURIER_MES1=clean_genes(cout[order(-mes1),head(gene_names,50)]),
  GAROFANO_GPM=clean_genes(garofano$GPM),
  NOMURA_MP_6_MES=clean_genes(nomura$MP_6_MES)
)
source_file <- c(
  CANONICAL_NEFTEL_MES_LIKE_95_ANCHOR=SIGREG,
  VERHAAK_MES=file.path(ICB,"GBM_Verhaak_signatures_df.rds"),
  NEFTEL_MES1=file.path(ICB,"Neftel_et_al_signatures.xlsx"),
  NEFTEL_MES2=file.path(ICB,"Neftel_et_al_signatures.xlsx"),
  COUTURIER_MES1=file.path(ICB,"Couturier_et_al_signatures.csv"),
  GAROFANO_GPM=file.path(ICB,"Garofano_et_al_signatures.xlsx"),
  NOMURA_MP_6_MES=file.path(ICB,"Nomura_et_al_signatures.xlsx")
)
source_pub <- c(
  CANONICAL_NEFTEL_MES_LIKE_95_ANCHOR="Frozen project canonical Neftel MES-like union",
  VERHAAK_MES="Verhaak et al.", NEFTEL_MES1="Neftel et al.", NEFTEL_MES2="Neftel et al.",
  COUTURIER_MES1="Couturier et al.", GAROFANO_GPM="Garofano et al.", NOMURA_MP_6_MES="Nomura et al."
)
universe <- Reduce(intersect,list(colnames(bulk),colnames(mal),colnames(eco)))
mapped <- lapply(sets,intersect,y=universe)
classifier_registry <- rbindlist(lapply(names(sets),function(id){
  n0 <- length(sets[[id]]); nm <- length(mapped[[id]])
  data.table(classifier=id,source_publication=source_pub[id],source_file=source_file[id],genes_original=n0,genes_mapped=nm,
             mapping_rate=nm/n0,score_definition="equal-weight sum of gene-wise z scores on log2(expression+1), divided by original certified gene count; same mapped genes in bulk/malignant/ecological",
             previously_used_in_project=TRUE,
             certification_status=if(nm/n0>=.8)"CERTIFIED_EVALUABLE" else "CLASSIFIER_NOT_EVALUABLE_MAPPING_BELOW_FROZEN_80_PERCENT_GATE",
             mixed_status=if(id=="CANONICAL_NEFTEL_MES_LIKE_95_ANCHOR")"FROZEN_STABLE_MIXED_AVAILABLE" else "SOURCE_DECOMPOSITION_NOT_EVALUABLE_FOR_CLASSIFIER")
}))
w(classifier_registry,"P1_CLASSIFIER_GENESET_REGISTRY.tsv")
eval_ids <- classifier_registry[certification_status=="CERTIFIED_EVALUABLE",classifier]

score_expr <- function(e, genes, denom) {
  x <- log2(e+1); z <- scale(x); z[!is.finite(z)] <- 0
  rowSums(z[,genes,drop=FALSE])/denom
}
scores <- copy(d0[,.(participant_id,sample_id,time_years,event,age_decade,setting_recurrent,stable_mixed_contribution_sd)])
score_long <- list()
for(id in eval_ids){
  g <- mapped[[id]]; den <- length(sets[[id]])
  b <- score_expr(bulk,g,den); m <- score_expr(mal,g,den); e <- score_expr(eco,g,den)
  scores[,(paste0(id,"__bulk")):=znum(b)]
  scores[,(paste0(id,"__malignant")):=znum(m)]
  scores[,(paste0(id,"__ecological")):=znum(e)]
  score_long[[id]] <- data.table(participant_id=d0$participant_id,classifier=id,bulk_score=b,malignant_score=m,ecological_score=e,
                                 bulk_score_sd=znum(b),malignant_score_sd=znum(m),ecological_score_sd=znum(e),
                                 mixed_score_sd=if(id=="CANONICAL_NEFTEL_MES_LIKE_95_ANCHOR")d0$stable_mixed_contribution_sd else NA_real_)
}
score_long <- rbindlist(score_long)
w(score_long,"P1_CLASSIFIER_PATIENT_SCORES.tsv")

boot_model <- function(form, dat, target_terms, B=1000L){
  f0 <- coxph(form,dat,x=TRUE,y=TRUE)
  app <- calc_c(dat,predict(f0,type="lp")); coefs <- matrix(NA_real_,B,length(target_terms),dimnames=list(NULL,target_terms)); opt <- rep(NA_real_,B)
  for(b in seq_len(B)){
    ii <- sample.int(nrow(dat),nrow(dat),TRUE); tr <- dat[ii]
    fb <- tryCatch(suppressWarnings(coxph(form,tr,x=TRUE,y=TRUE,control=coxph.control(iter.max=60))),error=function(e)NULL)
    if(is.null(fb)) next
    cf <- coef(fb); if(all(target_terms%in%names(cf))&&all(is.finite(cf[target_terms]))) coefs[b,] <- cf[target_terms]
    ct <- calc_c(tr,predict(fb,newdata=tr,type="lp")); co <- calc_c(dat,predict(fb,newdata=dat,type="lp")); if(is.finite(ct)&&is.finite(co)) opt[b] <- ct-co
  }
  slope <- tryCatch(unname(coef(coxph(Surv(time_years,event)~lp,data=data.frame(time_years=dat$time_years,event=dat$event,lp=predict(f0,type="lp"))))[1]),error=function(e)NA_real_)
  list(fit=f0,apparent_c=app,corrected_c=app-mean(opt,na.rm=TRUE),optimism=mean(opt,na.rm=TRUE),calibration_slope=slope,
       coef_boot=coefs,bootstrap_success=sum(complete.cases(coefs)),optimism_success=sum(is.finite(opt)))
}

set.seed(2026081201L)
model_rows <- list(); tv_rows <- list(); model_i <- 0L; tv_i <- 0L
for(id in eval_ids){
  q <- data.table(time_years=scores$time_years,event=scores$event,age_decade=scores$age_decade,setting_recurrent=scores$setting_recurrent,
                  bulk=scores[[paste0(id,"__bulk")]],malignant=scores[[paste0(id,"__malignant")]],ecological=scores[[paste0(id,"__ecological")]],
                  mixed=if(id=="CANONICAL_NEFTEL_MES_LIKE_95_ANCHOR")scores$stable_mixed_contribution_sd else NA_real_)
  specs <- list(
    M_bulk=list(form=Surv(time_years,event)~age_decade+setting_recurrent+bulk,terms="bulk"),
    M_malignant=list(form=Surv(time_years,event)~age_decade+setting_recurrent+malignant,terms="malignant"),
    M_ecological=list(form=Surv(time_years,event)~age_decade+setting_recurrent+ecological,terms="ecological"),
    M_joint=list(form=Surv(time_years,event)~age_decade+setting_recurrent+malignant+ecological,terms=c("malignant","ecological"))
  )
  if(id=="CANONICAL_NEFTEL_MES_LIKE_95_ANCHOR") specs$M_mixed <- list(form=Surv(time_years,event)~age_decade+setting_recurrent+mixed,terms="mixed")
  for(mid in names(specs)){
    z <- tryCatch(boot_model(specs[[mid]]$form,q,specs[[mid]]$terms,1000L),error=function(e)NULL)
    if(is.null(z)){
      model_i<-model_i+1L;model_rows[[model_i]]<-data.table(classifier=id,model_id=mid,term=NA_character_,status="NOT_ESTIMABLE")
      next
    }
    f <- z$fit; sm <- summary(f); ci <- exp(confint(f)); dfb <- tryCatch(max(abs(residuals(f,type="dfbeta")),na.rm=TRUE),error=function(e)NA_real_)
    for(te in specs[[mid]]$terms){
      bvec <- z$coef_boot[,te]; bvec <- bvec[is.finite(bvec)]
      model_i <- model_i+1L
      model_rows[[model_i]] <- data.table(classifier=id,model_id=mid,term=te,status="ESTIMABLE",n=f$n,events=f$nevent,beta=coef(f)[te],hr=exp(coef(f)[te]),
        ci_low=ci[te,1],ci_high=ci[te,2],p_value=sm$coefficients[te,"Pr(>|z|)"],ph_global_p=safe_ph(f),ph_term_p=safe_ph(f,te),
        apparent_c_index=z$apparent_c,optimism=z$optimism,optimism_corrected_c_index=z$corrected_c,calibration_slope=z$calibration_slope,
        bootstrap_success=length(bvec),bootstrap_favorable_fraction=mean(bvec<0),bootstrap_beta_median=median(bvec),bootstrap_beta_ci_low=quantile(bvec,.025),bootstrap_beta_ci_high=quantile(bvec,.975),max_abs_dfbeta=dfb)
      # Per protocol, any global or term violation triggers a fixed time-varying sensitivity for this component.
      if((is.finite(safe_ph(f))&&safe_ph(f)<.05)||(is.finite(safe_ph(f,te))&&safe_ph(f,te)<.05)){
        others <- setdiff(specs[[mid]]$terms,te)
        rhs <- paste(c("age_decade","setting_recurrent",others,te,paste0("tt(",te,")")),collapse="+")
        ft <- tryCatch(suppressWarnings(coxph(as.formula(paste0("Surv(time_years,event)~",rhs)),q,tt=function(x,t,...)x*log(pmax(t,1e-6)),x=TRUE,y=TRUE)),error=function(e)NULL)
        if(!is.null(ft)){
          cf<-coef(ft); vc<-vcov(ft); it<-paste0("tt(",te,")")
          for(tm in c(3,6,12,18)){
            L<-rep(0,length(cf));names(L)<-names(cf);L[te]<-1;L[it]<-log(tm/12);bb<-sum(L*cf);se<-sqrt(drop(t(L)%*%vc%*%L))
            tv_i<-tv_i+1L;tv_rows[[tv_i]]<-data.table(classifier=id,source_model=mid,component=te,time_months=tm,beta=bb,hr=exp(bb),ci_low=exp(bb-1.96*se),ci_high=exp(bb+1.96*se),
              trigger_global_ph_p=safe_ph(f),trigger_term_ph_p=safe_ph(f,te),interaction_beta=cf[it],interaction_p=summary(ft)$coefficients[it,"Pr(>|z|)"],status="TIME_VARYING_EVALUATED")
          }
        }
      } else {
        for(tm in c(3,6,12,18)){tv_i<-tv_i+1L;tv_rows[[tv_i]]<-data.table(classifier=id,source_model=mid,component=te,time_months=tm,beta=coef(f)[te],hr=exp(coef(f)[te]),ci_low=ci[te,1],ci_high=ci[te,2],trigger_global_ph_p=safe_ph(f),trigger_term_ph_p=safe_ph(f,te),interaction_beta=NA_real_,interaction_p=NA_real_,status="PH_COMPATIBLE_CONSTANT_EFFECT")}
      }
    }
  }
  if(id!="CANONICAL_NEFTEL_MES_LIKE_95_ANCHOR"){
    model_i<-model_i+1L;model_rows[[model_i]]<-data.table(classifier=id,model_id="M_mixed",term="mixed",status="SOURCE_DECOMPOSITION_NOT_EVALUABLE_FOR_CLASSIFIER")
  }
}
models <- rbindlist(model_rows,fill=TRUE); tv <- rbindlist(tv_rows,fill=TRUE)
w(models,"P1_CLASSIFIER_COMPONENT_MODELS.tsv"); w(tv,"P1_CLASSIFIER_TIME_VARYING_EFFECTS.tsv")

# Paired patient bootstrap: all classifiers re-fit on the same resample.
set.seed(2026081202L)
pb <- vector("list",1000L*length(eval_ids)); kk <- 0L
for(b in seq_len(1000L)){
  ii <- sample.int(nrow(scores),nrow(scores),TRUE)
  for(id in eval_ids){
    q <- data.table(time_years=scores$time_years[ii],event=scores$event[ii],age_decade=scores$age_decade[ii],setting_recurrent=scores$setting_recurrent[ii],
                    malignant=scores[[paste0(id,"__malignant")]][ii],ecological=scores[[paste0(id,"__ecological")]][ii],
                    mixed=if(id=="CANONICAL_NEFTEL_MES_LIKE_95_ANCHOR")scores$stable_mixed_contribution_sd[ii] else NA_real_)
    fm <- tryCatch(suppressWarnings(coxph(Surv(time_years,event)~age_decade+setting_recurrent+malignant,q)),error=function(e)NULL)
    fe <- tryCatch(suppressWarnings(coxph(Surv(time_years,event)~age_decade+setting_recurrent+ecological,q)),error=function(e)NULL)
    fx <- if(id=="CANONICAL_NEFTEL_MES_LIKE_95_ANCHOR")tryCatch(suppressWarnings(coxph(Surv(time_years,event)~age_decade+setting_recurrent+mixed,q)),error=function(e)NULL) else NULL
    bm<-safe_coef(fm,"malignant");be<-safe_coef(fe,"ecological");bx<-safe_coef(fx,"mixed")
    kk<-kk+1L; pb[[kk]]<-data.table(iteration=b,classifier=id,beta_malignant=bm,beta_ecological=be,beta_mixed=bx,
      ecological_more_favorable_than_malignant=is.finite(be)&is.finite(bm)&be<bm,ecological_favorable=is.finite(be)&be<0,malignant_favorable=is.finite(bm)&bm<0,
      mixed_favorable=if(is.finite(bx))bx<0 else NA,models_estimable=is.finite(be)&is.finite(bm),
      ph_compatible=is.finite(safe_ph(fm))&is.finite(safe_ph(fe))&safe_ph(fm)>=.05&safe_ph(fe)>=.05)
  }
}
pb <- rbindlist(pb[seq_len(kk)]); w(pb,"P1_CLASSIFIER_SOURCE_DOMINANCE_BOOTSTRAP.tsv")

cons <- pb[,.(bootstrap_attempts=.N,joint_estimable=sum(models_estimable),ecological_favorable_fraction=mean(ecological_favorable[models_estimable]),
               malignant_favorable_fraction=mean(malignant_favorable[models_estimable]),ecological_more_favorable_fraction=mean(ecological_more_favorable_than_malignant[models_estimable]),
               mixed_favorable_fraction=if(all(is.na(mixed_favorable)))NA_real_ else mean(mixed_favorable,na.rm=TRUE),ph_compatible_fraction=mean(ph_compatible[models_estimable])),by=classifier]
for(tm in c(6,12,18)){
  et <- tv[source_model=="M_ecological" & component=="ecological" & time_months==tm,.(classifier,eco_beta=beta)]
  mt <- tv[source_model=="M_malignant" & component=="malignant" & time_months==tm,.(classifier,mal_beta=beta)]
  z <- merge(et,mt,by="classifier",all=TRUE); z[,direction:=fifelse(eco_beta<0&eco_beta<mal_beta,"ECOLOGICAL_FAVORABLE_AND_MORE_FAVORABLE",fifelse(eco_beta<0,"ECOLOGICAL_FAVORABLE_NOT_MORE_FAVORABLE","ECOLOGICAL_NOT_FAVORABLE"))]
  setnames(z,c("eco_beta","mal_beta","direction"),paste0(c("ecological_beta_","malignant_beta_","direction_"),tm,"m")); cons<-merge(cons,z,by="classifier",all.x=TRUE)
}
cons[,criteria_ecological_favorable_ge80:=ecological_favorable_fraction>=.8]
cons[,criteria_ecological_more_favorable_ge80:=ecological_more_favorable_fraction>=.8]
cons[,criteria_no_systematic_reversal_6_12_18:=direction_6m!="ECOLOGICAL_NOT_FAVORABLE"&direction_12m!="ECOLOGICAL_NOT_FAVORABLE"&direction_18m!="ECOLOGICAL_NOT_FAVORABLE"]
cons[,classifier_gate_pass:=criteria_ecological_favorable_ge80&criteria_ecological_more_favorable_ge80&criteria_no_systematic_reversal_6_12_18]
summary_cons <- data.table(classifier="ACROSS_CLASSIFIER_SUMMARY",bootstrap_attempts=NA_integer_,joint_estimable=sum(cons$classifier_gate_pass),
  ecological_favorable_fraction=mean(cons$criteria_ecological_favorable_ge80),malignant_favorable_fraction=NA_real_,ecological_more_favorable_fraction=mean(cons$criteria_ecological_more_favorable_ge80),
  mixed_favorable_fraction=NA_real_,ph_compatible_fraction=NA_real_,criteria_ecological_favorable_ge80=sum(cons$criteria_ecological_favorable_ge80),
  criteria_ecological_more_favorable_ge80=sum(cons$criteria_ecological_more_favorable_ge80),criteria_no_systematic_reversal_6_12_18=sum(cons$criteria_no_systematic_reversal_6_12_18),
  classifier_gate_pass=sum(cons$classifier_gate_pass))
w(rbindlist(list(cons,summary_cons),fill=TRUE),"P1_CROSS_CLASSIFIER_CONSISTENCY.tsv")

# P2 benchmark registry and values.
unified <- fread(file.path(ICB,"bulk_transcriptional_classifiers/icb_cohort_unified.tsv"))
sc <- fread(file.path(ICB,"bulk_transcriptional_classifiers/sc_classifier_scores.tsv"))
bm <- copy(d0[,.(participant_id,sample_id,time_years,event,age_decade,setting_recurrent,stable_mixed_contribution_sd,ecological_component_sd,malignant_fraction,myeloid_fraction,vascular_stromal_fraction)])
bm[,`:=`(
  Neftel_T_cell_score=sc$Neftel_t_cell[match(sample_id,sc$sample_id)],
  Neftel_macrophage_score=sc$Neftel_macrophage[match(sample_id,sc$sample_id)],
  HLA_I_ssGSEA=unified$GOCC_MHC_CLASS_I_PROTEIN_COMPLEX_pre[match(participant_id,unified$participant_id)],
  HLA_II_ssGSEA=unified$GOCC_MHC_CLASS_II_PROTEIN_COMPLEX_pre[match(participant_id,unified$participant_id)]
)]
getsum <- function(pattern){cc<-grep(pattern,colnames(theta),ignore.case=TRUE,value=TRUE);if(length(cc))rowSums(theta[,cc,drop=FALSE]) else rep(NA_real_,nrow(theta))}
bm[,T_cell_fraction:=getsum("t.?cell|lymph")[sample_id]]

readgmt <- function(id){ln<-readLines(GMT); hit<-ln[startsWith(ln,paste0(id,"\t"))]; if(length(hit)!=1)stop("GMT ID not unique: ",id); toupper(strsplit(hit,"\t")[[1]][-c(1,2)])}
bench_sets <- list(
  IFN_gamma_hallmark=readgmt("HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  hypoxia_hallmark=unique(reg[signature_id=="HALLMARK_HYPOXIA",gene]),
  matrix_core_matrisome=unique(reg[signature_id=="NABA_CORE_MATRISOME",gene])
)
bulk_log <- log2(bulk+1); bulk_z <- scale(bulk_log); bulk_z[!is.finite(bulk_z)]<-0
for(nm in names(bench_sets)){g<-intersect(bench_sets[[nm]],colnames(bulk_z));bm[,(nm):=rowMeans(bulk_z[,g,drop=FALSE])]}
pathology <- fread(file.path(FORMAL,"04_orthogonal/HLA_CD3_PATIENT_LEVEL_TABLE.tsv"))
path_cols <- intersect(c("participant_id","cd3_density","tumor_cell_HLA_I","tumor_cell_HLA_II"),names(pathology))
bm <- merge(bm,pathology[,..path_cols],by="participant_id",all.x=TRUE)
setnames(bm,old=intersect(c("cd3_density","tumor_cell_HLA_I","tumor_cell_HLA_II"),names(bm)),new=c("pathology_CD3_density","pathology_tumor_HLA_I","pathology_tumor_HLA_II")[match(intersect(c("cd3_density","tumor_cell_HLA_I","tumor_cell_HLA_II"),names(bm)),c("cd3_density","tumor_cell_HLA_I","tumor_cell_HLA_II"))])

bench_info <- data.table(
  benchmark=c("Neftel_T_cell_score","IFN_gamma_hallmark","HLA_I_ssGSEA","HLA_II_ssGSEA","Neftel_macrophage_score","T_cell_fraction","myeloid_fraction","malignant_fraction","vascular_stromal_fraction","hypoxia_hallmark","matrix_core_matrisome","pathology_CD3_density","pathology_tumor_HLA_I","pathology_tumor_HLA_II"),
  source=c("Official ICB sc_classifier_scores.tsv","Frozen V1.5 MSigDB Hallmark","Official ICB GOCC ssGSEA","Official ICB GOCC ssGSEA","Official ICB sc_classifier_scores.tsv","Frozen V1.5 BayesPrism theta","Frozen V1.5 BayesPrism theta","Frozen V1.5 BayesPrism theta","Frozen V1.5 BayesPrism theta","Frozen project signature registry","Frozen project signature registry","R7 pathology table","R7 pathology table","R7 pathology table"),
  type=c(rep("continuous",14)),
  score_definition=c("Official project score","mean gene-z; frozen Hallmark genes","Official ssGSEA","Official ssGSEA","Official project score","posterior fraction","posterior fraction","posterior fraction / purity proxy","posterior fraction","mean gene-z","mean gene-z","cells/mm2","tumor-cell HLA-I readout","tumor-cell HLA-II readout"),
  certified=TRUE,
  priority=c("PRIMARY","PRIMARY","SECONDARY","PRIMARY","SECONDARY","SECONDARY","PRIMARY","PRIMARY","SECONDARY","CONTEXT_CONTROL","CONTEXT_CONTROL","ORTHOGONAL_LOW_N","ORTHOGONAL_LOW_N","ORTHOGONAL_LOW_N")
)
bench_info[,denominator:=vapply(benchmark,function(v)sum(is.finite(bm[[v]])),integer(1))]
bench_info[,missingness:=40-denominator]
w(bench_info,"P2_BENCHMARK_AVAILABILITY.tsv")
w(bm,"P2_BENCHMARK_PATIENT_TABLE.tsv")

fit_att <- function(q, score, benchmark, B=1000L, seed=1L){
  set.seed(seed)
  q <- copy(q); q[,score:=get(score)]; q[,benchmark:=znum(get(benchmark))]
  q <- q[complete.cases(q[,.(time_years,event,age_decade,setting_recurrent,score,benchmark)])]
  f0<-coxph(Surv(time_years,event)~age_decade+setting_recurrent+score,q,x=TRUE,y=TRUE)
  f1<-coxph(Surv(time_years,event)~age_decade+setting_recurrent+score+benchmark,q,x=TRUE,y=TRUE)
  cf0<-coef(f0)["score"];cf1<-coef(f1)["score"];ci0<-exp(confint(f0)["score",]);ci1<-exp(confint(f1)["score",])
  boot0<-boot1<-rep(NA_real_,B)
  for(b in seq_len(B)){ii<-sample.int(nrow(q),nrow(q),TRUE);x<-q[ii];g0<-tryCatch(suppressWarnings(coxph(Surv(time_years,event)~age_decade+setting_recurrent+score,x)),error=function(e)NULL);g1<-tryCatch(suppressWarnings(coxph(Surv(time_years,event)~age_decade+setting_recurrent+score+benchmark,x)),error=function(e)NULL);boot0[b]<-safe_coef(g0,"score");boot1[b]<-safe_coef(g1,"score")}
  r<-suppressWarnings(cor(q$score,q$benchmark));vif<-1/(1-r^2);cond<-kappa(cbind(1,q$age_decade,q$setting_recurrent,q$score,q$benchmark),exact=TRUE)
  sm1<-summary(f1)
  data.table(n=nrow(q),events=sum(q$event),beta_base=cf0,hr_base=exp(cf0),ci_base_low=ci0[1],ci_base_high=ci0[2],p_base=summary(f0)$coefficients["score","Pr(>|z|)"],
    beta_adjusted=cf1,hr_adjusted=exp(cf1),ci_adjusted_low=ci1[1],ci_adjusted_high=ci1[2],p_adjusted=sm1$coefficients["score","Pr(>|z|)"],
    attenuation_signed=1-cf1/cf0,attenuation_percent=100*(1-cf1/cf0),sign_preserved=sign(cf0)==sign(cf1),bootstrap_base_success=sum(is.finite(boot0)),bootstrap_adjusted_success=sum(is.finite(boot1)),
    bootstrap_base_favorable=mean(boot0<0,na.rm=TRUE),bootstrap_adjusted_favorable=mean(boot1<0,na.rm=TRUE),benchmark_beta=coef(f1)["benchmark"],benchmark_p=sm1$coefficients["benchmark","Pr(>|z|)"],
    correlation=r,VIF=vif,condition_number=cond,separability_status=if(vif>5||cond>30)"NON_SEPARABLE_DUE_TO_COLLINEARITY" else "SEPARABLE",
    ph_base_global_p=safe_ph(f0),ph_base_score_p=safe_ph(f0,"score"),ph_adjusted_global_p=safe_ph(f1),ph_adjusted_score_p=safe_ph(f1,"score"))
}

set.seed(2026081203L)
stable_rows <- list(); eco_rows <- list(); ai<-0L;ei<-0L
for(j in seq_len(nrow(bench_info))){
  v<-bench_info$benchmark[j];a<-fit_att(bm,"stable_mixed_contribution_sd",v,1000L,2026081203L+j);a[,benchmark:=v];a[,priority:=bench_info$priority[j]];ai<-ai+1L;stable_rows[[ai]]<-a
  # Ecological constant model first; if score/global PH violates, use fixed-time TV coefficients.
  q<-copy(bm);q[,score:=ecological_component_sd];q[,benchmark:=znum(get(v))];q<-q[complete.cases(q[,.(time_years,event,age_decade,setting_recurrent,score,benchmark)])]
  f0<-coxph(Surv(time_years,event)~age_decade+setting_recurrent+score,q,x=TRUE,y=TRUE);f1<-coxph(Surv(time_years,event)~age_decade+setting_recurrent+score+benchmark,q,x=TRUE,y=TRUE)
  trigger<-(is.finite(safe_ph(f0))&&safe_ph(f0)<.05)||(is.finite(safe_ph(f0,"score"))&&safe_ph(f0,"score")<.05)||(is.finite(safe_ph(f1))&&safe_ph(f1)<.05)||(is.finite(safe_ph(f1,"score"))&&safe_ph(f1,"score")<.05)
  times<-if(trigger)c(6,12,18) else NA_real_
  if(trigger){
    ft0<-coxph(Surv(time_years,event)~age_decade+setting_recurrent+score+tt(score),q,tt=function(x,t,...)x*log(pmax(t,1e-6)),x=TRUE,y=TRUE)
    ft1<-coxph(Surv(time_years,event)~age_decade+setting_recurrent+score+benchmark+tt(score),q,tt=function(x,t,...)x*log(pmax(t,1e-6)),x=TRUE,y=TRUE)
    for(tm in times){
      beta_t<-function(f){cf<-coef(f);cf["score"]+cf["tt(score)"]*log(tm/12)}
      b0<-beta_t(ft0);b1<-beta_t(ft1);bb0<-bb1<-rep(NA_real_,1000L)
      set.seed(2026082200L+j*10+tm)
      for(b in seq_len(1000L)){ii<-sample.int(nrow(q),nrow(q),TRUE);x<-q[ii];g0<-tryCatch(suppressWarnings(coxph(Surv(time_years,event)~age_decade+setting_recurrent+score+tt(score),x,tt=function(z,t,...)z*log(pmax(t,1e-6)))),error=function(e)NULL);g1<-tryCatch(suppressWarnings(coxph(Surv(time_years,event)~age_decade+setting_recurrent+score+benchmark+tt(score),x,tt=function(z,t,...)z*log(pmax(t,1e-6)))),error=function(e)NULL);if(!is.null(g0))bb0[b]<-beta_t(g0);if(!is.null(g1))bb1[b]<-beta_t(g1)}
      r<-cor(q$score,q$benchmark);vif<-1/(1-r^2);cond<-kappa(cbind(1,q$age_decade,q$setting_recurrent,q$score,q$benchmark),exact=TRUE)
      ei<-ei+1L;eco_rows[[ei]]<-data.table(benchmark=v,priority=bench_info$priority[j],n=nrow(q),events=sum(q$event),effect_summary="TIME_VARYING_BETA",time_months=tm,beta_base=b0,beta_adjusted=b1,hr_base=exp(b0),hr_adjusted=exp(b1),attenuation_signed=1-b1/b0,attenuation_percent=100*(1-b1/b0),sign_preserved=sign(b0)==sign(b1),bootstrap_base_success=sum(is.finite(bb0)),bootstrap_adjusted_success=sum(is.finite(bb1)),bootstrap_base_favorable=mean(bb0<0,na.rm=TRUE),bootstrap_adjusted_favorable=mean(bb1<0,na.rm=TRUE),correlation=r,VIF=vif,condition_number=cond,separability_status=if(vif>5||cond>30)"NON_SEPARABLE_DUE_TO_COLLINEARITY" else "SEPARABLE",ph_base_global_p=safe_ph(f0),ph_base_score_p=safe_ph(f0,"score"),ph_adjusted_global_p=safe_ph(f1),ph_adjusted_score_p=safe_ph(f1,"score"))
    }
  } else {
    a<-fit_att(bm,"ecological_component_sd",v,1000L,2026083200L+j);a[,`:=`(benchmark=v,priority=bench_info$priority[j],effect_summary="CONSTANT_BETA",time_months=NA_real_)];ei<-ei+1L;eco_rows[[ei]]<-a
  }
}
stable_att<-rbindlist(stable_rows,fill=TRUE);stable_att[,benchmark_BH_FDR:=p.adjust(benchmark_p,method="BH")]
eco_att<-rbindlist(eco_rows,fill=TRUE)
w(stable_att,"P2_STABLE_MIXED_BENCHMARK_ATTENUATION.tsv");w(eco_att,"P2_ECOLOGICAL_BENCHMARK_ATTENUATION.tsv")

primary_bench <- c("Neftel_T_cell_score","IFN_gamma_hallmark","HLA_II_ssGSEA","myeloid_fraction","malignant_fraction")
pa <- stable_att[benchmark%in%primary_bench]
adjudication <- if(all(pa$sign_preserved)&mean(abs(pa$attenuation_percent)<30)>=.6&all(pa$bootstrap_adjusted_favorable>=.8)&all(pa$separability_status=="SEPARABLE")) {
  "STABLE_MIXED_NOT_FULLY_EXPLAINED_BY_GENERIC_IMMUNE_HOTNESS"
} else if(mean(pa$attenuation_percent>=30|!pa$sign_preserved)>=.5) {
  "STABLE_MIXED_CLINICAL_SIGNAL_LARGELY_EMBEDS_IMMUNE_ECOLOGICAL_CONTEXT"
} else {
  "IMMUNE_SPECIFICITY_NOT_RESOLVABLE_DUE_TO_COLLINEARITY_OR_DENOMINATOR"
}
adj <- data.table(adjudication=adjudication,primary_benchmark_n=nrow(pa),favorable_direction_retained=sum(pa$sign_preserved),attenuation_below_30_percent=sum(abs(pa$attenuation_percent)<30),bootstrap_favorable_ge80=sum(pa$bootstrap_adjusted_favorable>=.8),separable=sum(pa$separability_status=="SEPARABLE"),
                  largest_attenuation_benchmark=stable_att[which.max(attenuation_percent),benchmark],largest_attenuation_percent=max(stable_att$attenuation_percent,na.rm=TRUE),
                  largest_primary_attenuation_benchmark=pa[which.max(attenuation_percent),benchmark],largest_primary_attenuation_percent=max(pa$attenuation_percent,na.rm=TRUE),
                  claim_boundary="Association specificity within a treated cohort; not independent treatment prediction")
w(adj,"P2_IMMUNE_SPECIFICITY_ADJUDICATION.tsv")

writeLines(capture.output(sessionInfo()),file.path(RUN,"logs/PRIMARY_R_SESSION_INFO.txt"))
cat("PRIMARY_SPECIFICITY_COMPLETE",length(eval_ids),nrow(models),nrow(pb),adjudication,"\n")
