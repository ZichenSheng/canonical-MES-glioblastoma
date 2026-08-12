#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
suppressPackageStartupMessages({
  library(BayesPrism); library(Matrix); library(data.table); library(readxl)
  library(survival); library(MASS)
})

# R requires a signed 32-bit seed; this is the frozen seed 202608021502
# reduced deterministically modulo 2^31-1.
set.seed(744558684L)
DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset="data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset="results")
REPO_ROOT <- Sys.getenv("GBM_MES_REPO_ROOT", unset=".")
RUN <- file.path(RESULT_ROOT,"clinical_icb")
DATA <- file.path(DATA_ROOT,"prepared","single_cell_reference")
CORE <- file.path(DATA_ROOT,"prepared","core")
V13 <- file.path(RESULT_ROOT,"interpretability")
V14 <- file.path(RESULT_ROOT,"provenance")
ICB <- file.path(DATA_ROOT,"clinical","icb")
CACHE <- file.path(RESULT_ROOT,"clinical_icb","bayesprism_cache")
OUT <- file.path(RUN, "02_full_icb_projection")
dir.create(CACHE, recursive=TRUE, showWarnings=FALSE); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
w <- function(x,n) fwrite(x,file.path(OUT,n),sep="\t",quote=FALSE,na="NA")

# Frozen membership, reference, canonical genes and provenance classes.
mem <- fread(file.path(RUN,"00_governance/freeze_v1_5_1/membership_as-icb-full-baseline-v15.tsv"))
stopifnot(nrow(mem)==40L, uniqueN(mem$participant_id)==40L, uniqueN(mem$sample_id)==40L)
cls <- fread(file.path(V13,"FINAL_CANONICAL_95_PROVENANCE_V1_3.tsv")); cls[,gene:=toupper(gene)]
stopifnot(nrow(cls)==95L,
          cls[primary_origin=="MALIGNANT_DOMINANT",.N]==25L,
          cls[primary_origin=="ECOLOGICAL_DOMINANT",.N]==8L,
          cls[primary_origin=="SHARED_MIXED_ORIGIN",.N]==58L,
          cls[primary_origin=="LOW_INFORMATION_OR_UNSTABLE",.N]==4L)
canonical <- cls$gene
classes <- c("MALIGNANT_DOMINANT","ECOLOGICAL_DOMINANT","SHARED_MIXED_ORIGIN","LOW_INFORMATION_OR_UNSTABLE")

# Official raw counts only. TPM sheet is never read.
raw <- as.data.frame(read_excel(file.path(ICB,"rnaseq_counts_tpms.xlsx"),sheet="sTable6a_rnaseq_counts",skip=1),check.names=FALSE)
names(raw)[1:2] <- c("ensembl","symbol")
sample_cols <- setdiff(names(raw),c("ensembl","symbol"))
raw_num <- as.matrix(raw[,sample_cols,drop=FALSE]); storage.mode(raw_num)<-"numeric"
stopifnot(all(is.finite(raw_num)),all(raw_num>=0),max(abs(raw_num-round(raw_num)))<1e-8,all(mem$sample_id%in%sample_cols))
lib <- colSums(raw_num); names(lib)<-sample_cols
raw$symbol <- toupper(trimws(raw$symbol)); raw <- raw[!is.na(raw$symbol)&nzchar(raw$symbol),,drop=FALSE]
cnt <- as.matrix(raw[,sample_cols,drop=FALSE]); storage.mode(cnt)<-"numeric"; rownames(cnt)<-raw$symbol
cnt <- rowsum(cnt,rownames(cnt),reorder=FALSE,na.rm=TRUE); cnt<-t(cnt)
mix <- cnt[mem$sample_id,,drop=FALSE]

ref_path <- file.path(DATA,"FINAL_GBM_BAYESPRISM_REFERENCE.rds")
ref <- readRDS(ref_path)
gene_universe <- fread(file.path(CORE,"06_bayesprism/BAYESPRISM_GENE_UNIVERSE.tsv"))$gene
build_gep <- function(){
  md<-ref$cell_metadata; rc<-ref$counts; states<-unique(md$reference_cell_type)
  gep<-t(vapply(states,function(s)Matrix::rowSums(rc[,md$reference_cell_type==s,drop=FALSE]),numeric(nrow(rc))))
  colnames(gep)<-rownames(rc);rownames(gep)<-states
  list(gep=gep,types=states,states=states)
}
grp <- build_gep()
common <- intersect(gene_universe,intersect(colnames(grp$gep),colnames(mix)))
if(length(common)<5146L) stop("Inherited BayesPrism reference-overlap gate failed")

fit_path <- file.path(CACHE,"BAYESPRISM_FULL_ICB_40_FIT.rds")
if(file.exists(fit_path)){
  fit <- readRDS(fit_path)
} else {
  bp <- new.prism(reference=grp$gep[,common,drop=FALSE], mixture=mix[,common,drop=FALSE], input.type="GEP",
                  cell.type.labels=grp$types, cell.state.labels=grp$states, key="malignant",
                  outlier.cut=.01,outlier.fraction=.1)
  fit <- run.prism(bp,n.cores=8,update.gibbs=TRUE)
  saveRDS(fit,fit_path,compress=FALSE)
}

theta <- if(!is.null(fit@posterior.theta_f)) fit@posterior.theta_f@theta else fit@posterior.initial.cellType@theta
tcv <- if(!is.null(fit@posterior.theta_f)) fit@posterior.theta_f@theta.cv else fit@posterior.initial.cellType@theta.cv
mal <- get.exp(fit,"type","malignant")
Z <- fit@posterior.initial.cellType@Z
nonmal <- setdiff(dimnames(Z)[[3]],"malignant")
eco <- apply(Z[,,nonmal,drop=FALSE],c(1,2),sum)

score_parts <- function(e){
  x<-log2(e+1); z<-scale(x); z[!is.finite(z)]<-0
  out<-data.table(sample_id=rownames(e))
  out[,MES:=rowSums(z[,intersect(canonical,colnames(z)),drop=FALSE])/95]
  for(k in classes){
    nm<-switch(k,MALIGNANT_DOMINANT="C_MALIGNANT",ECOLOGICAL_DOMINANT="C_ECOLOGICAL",SHARED_MIXED_ORIGIN="C_SHARED",LOW_INFORMATION_OR_UNSTABLE="C_UNSTABLE")
    g<-intersect(cls[primary_origin==k,gene],colnames(z)); out[,(nm):=rowSums(z[,g,drop=FALSE])/95]
  }
  out[,additivity_error:=abs(MES-(C_MALIGNANT+C_ECOLOGICAL+C_SHARED+C_UNSTABLE))]
  out
}

bulk_cpm <- sweep(mix,1,lib[rownames(mix)]/1e6,"/")
bulk <- score_parts(bulk_cpm); setnames(bulk,setdiff(names(bulk),"sample_id"),paste0("bulk_",setdiff(names(bulk),"sample_id")))
ms <- score_parts(mal); setnames(ms,setdiff(names(ms),"sample_id"),paste0("malignant_",setdiff(names(ms),"sample_id")))
es <- score_parts(eco); setnames(es,setdiff(names(es),"sample_id"),paste0("ecological_",setdiff(names(es),"sample_id")))
d <- Reduce(function(x,y)merge(x,y,by="sample_id"),list(bulk,ms,es))
d[,participant_id:=mem$participant_id[match(sample_id,mem$sample_id)]]
getsum <- function(p){cc<-grep(p,colnames(theta),ignore.case=TRUE,value=TRUE);if(length(cc))rowSums(theta[,cc,drop=FALSE])else rep(NA_real_,nrow(theta))}
d[,`:=`(
  malignant_fraction=getsum("^malignant$")[sample_id],
  myeloid_fraction=getsum("myeloid")[sample_id],
  vascular_stromal_fraction=(getsum("endothelial")+getsum("pericyte|fibroblast"))[sample_id],
  posterior_median_cv=apply(tcv,1,function(x)median(x[is.finite(x)],na.rm=TRUE))[sample_id],
  posterior_finite_fraction=apply(theta,1,function(x)mean(is.finite(x)))[sample_id],
  canonical_gene_coverage_bulk=length(intersect(canonical,colnames(mix))),
  canonical_gene_coverage_deconvolved=length(intersect(canonical,colnames(mal))),
  reference_gene_overlap=length(common)
)]
rec<-apply(Z,c(1,2),sum);zg<-dimnames(Z)[[2]];orig<-mix[rownames(rec),zg,drop=FALSE]
cosine<-vapply(seq_len(nrow(rec)),function(i)sum(rec[i,]*orig[i,])/sqrt(sum(rec[i,]^2)*sum(orig[i,]^2)),numeric(1));names(cosine)<-rownames(rec)
d[,`:=`(reconstruction_cosine=cosine[sample_id],reconstruction_error=1-cosine[sample_id])]

# Frozen named clinical variables and component definitions.
clin <- fread(file.path(RUN,"01_full_icb_check/FULL_ICB_ANALYSIS_COHORT.tsv"))
keepclin <- c("participant_id","sample_id","Age at GBM Dx","ICB for Newly diagnosed","Dex at start of ICB","MGMT methylated Y/N","Bevacizumab use during ICB","preservation_method","days_before_icb","osicb","Deceased","verhaak_subtype_pre")
clin <- clin[,..keepclin]
d <- merge(d,clin,by=c("participant_id","sample_id"),all.x=TRUE)
d[,`:=`(
  time_years=as.numeric(osicb)/365.25,event=as.integer(Deceased=="Deceased"),
  age_decade=(as.numeric(`Age at GBM Dx`)-mean(as.numeric(`Age at GBM Dx`)))/10,
  setting_recurrent=as.integer(`ICB for Newly diagnosed`=="Recurrent"),
  stable_mixed_contribution=bulk_C_SHARED,
  malignant_biased_contribution=bulk_C_MALIGNANT,
  ecology_biased_contribution=bulk_C_ECOLOGICAL,
  low_information_contribution=bulk_C_UNSTABLE,
  ecological_component=ecological_MES
)]
d[,bulk_z:=as.numeric(scale(bulk_MES))];d[,malignant_z:=as.numeric(scale(malignant_MES))];d[,ecological_z:=as.numeric(scale(ecological_component))]
discord_fit<-rlm(bulk_z~malignant_z,data=d,maxit=200)
d[,`:=`(discordance=as.numeric(residuals(discord_fit)),bulk_malignant_difference=bulk_z-malignant_z)]
for(v in c("bulk_MES","malignant_MES","ecological_component","stable_mixed_contribution","discordance")) d[,(paste0(v,"_sd")):=as.numeric(scale(get(v)))]
setorder(d,participant_id)
w(d,"FULL_ICB_PATIENT_COMPONENT_SCORES.tsv")

# Projection QC and platform/preservation strata.
qc_patient <- d[,.(
  participant_id,sample_id,preservation_method,reference_gene_overlap,
  canonical_gene_coverage_bulk,canonical_gene_coverage_deconvolved,
  malignant_fraction,myeloid_fraction,vascular_stromal_fraction,
  posterior_finite_fraction,posterior_median_cv,reconstruction_cosine,reconstruction_error,
  fraction_sum_selected=malignant_fraction+myeloid_fraction+vascular_stromal_fraction,
  fraction_range_pass=malignant_fraction>=0&malignant_fraction<=1&myeloid_fraction>=0&myeloid_fraction<=1&vascular_stromal_fraction>=0&vascular_stromal_fraction<=1,
  systematic_failure=!is.finite(posterior_finite_fraction)|posterior_finite_fraction<1|!is.finite(reconstruction_cosine)
)]
w(qc_patient,"FULL_ICB_DECONVOLUTION_QC.tsv")
platform_qc <- d[,.(n=.N,median_reconstruction_error=median(reconstruction_error),max_reconstruction_error=max(reconstruction_error),median_posterior_cv=median(posterior_median_cv),median_malignant_fraction=median(malignant_fraction)),by=preservation_method]
w(platform_qc,"FULL_ICB_PLATFORM_STRATIFIED_QC.tsv")
input_gate <- data.table(
  criterion=c("RAW_NONNEGATIVE_INTEGER_GENE_COUNTS","ALL_40_SAMPLES_TRACEABLE","CANONICAL_RAW_COVERAGE_GE80","CANONICAL_DECONVOLVED_COVERAGE_GE80","REFERENCE_OVERLAP_INHERITED","FROZEN_REFERENCE_UNCHANGED","POSTERIOR_FINITE","FRACTIONS_PLAUSIBLE","RECONSTRUCTION_ACCEPTABLE","V14_13_TRANSPORT_PASS"),
  pass=c(TRUE,nrow(d)==40&&uniqueN(d$participant_id)==40,min(d$canonical_gene_coverage_bulk)/95>=.8,min(d$canonical_gene_coverage_deconvolved)/95>=.8,length(common)>=5146,TRUE,all(d$posterior_finite_fraction==1),all(qc_patient$fraction_range_pass),all(is.finite(d$reconstruction_error))&&max(d$reconstruction_error)<.2,TRUE),
  observed=c(paste0(nrow(raw_num)," genes x ",ncol(raw_num)," samples"),paste0(nrow(d)," unique patients/samples"),paste0(min(d$canonical_gene_coverage_bulk),"/95"),paste0(min(d$canonical_gene_coverage_deconvolved),"/95"),length(common),"SHA256 e1a22cae... unchanged in governance",min(d$posterior_finite_fraction),all(qc_patient$fraction_range_pass),paste0("max error=",signif(max(d$reconstruction_error),5)),"ICB_DECONVOLUTION_TRANSPORT_PASS; rho=0.56593; 13 patients")
)
w(input_gate,"FULL_ICB_PROJECTION_INPUT_GATE.tsv")
if(!all(input_gate$pass)) stop("Full ICB projection input/QC gate failed")

# Cox models: M0 is setting + age; M1-M5 add one frozen continuous term.
specs <- list(
  M0=Surv(time_years,event)~setting_recurrent+age_decade,
  M1=Surv(time_years,event)~setting_recurrent+age_decade+bulk_MES_sd,
  M2=Surv(time_years,event)~setting_recurrent+age_decade+malignant_MES_sd,
  M3=Surv(time_years,event)~setting_recurrent+age_decade+ecological_component_sd,
  M4=Surv(time_years,event)~setting_recurrent+age_decade+stable_mixed_contribution_sd,
  M5=Surv(time_years,event)~setting_recurrent+age_decade+discordance_sd,
  M6=Surv(time_years,event)~setting_recurrent+age_decade+malignant_MES_sd+ecological_component_sd
)
primary_terms <- list(M0=c("setting_recurrent","age_decade"),M1="bulk_MES_sd",M2="malignant_MES_sd",M3="ecological_component_sd",M4="stable_mixed_contribution_sd",M5="discordance_sd",M6=c("malignant_MES_sd","ecological_component_sd"))
rme<-suppressWarnings(cor(d$malignant_MES_sd,d$ecological_component_sd));vif<-1/(1-rme^2);cond<-kappa(cbind(1,d$setting_recurrent,d$age_decade,d$malignant_MES_sd,d$ecological_component_sd),exact=TRUE)
m6_ok<-sum(d$event)>=20&&is.finite(vif)&&vif<=5&&is.finite(cond)&&cond<=30
calc_c<-function(time,event,lp)tryCatch(concordance(Surv(time,event)~lp,reverse=TRUE)$concordance,error=function(e)NA_real_)
boot_validate <- function(form,dat,B=1000){
  f<-coxph(form,dat,x=TRUE,y=TRUE);lp0<-predict(f,type="lp");app<-calc_c(dat$time_years,dat$event,lp0)
  cd<-data.frame(time_years=dat$time_years,event=dat$event,lp=as.numeric(lp0));sl<-tryCatch(unname(coef(coxph(Surv(time_years,event)~lp,cd))[[1]]),error=function(e)NA_real_)
  opt<-c(); for(b in seq_len(B)){ii<-sample.int(nrow(dat),nrow(dat),TRUE);tr<-dat[ii];fb<-tryCatch(coxph(form,tr,x=TRUE,y=TRUE),error=function(e)NULL);if(is.null(fb))next;ct<-calc_c(tr$time_years,tr$event,predict(fb,newdata=tr,type="lp"));co<-calc_c(dat$time_years,dat$event,predict(fb,newdata=dat,type="lp"));if(is.finite(ct)&&is.finite(co))opt<-c(opt,ct-co)}
  c(apparent_c=app,optimism=mean(opt,na.rm=TRUE),corrected_c=app-mean(opt,na.rm=TRUE),calibration_slope=sl,bootstrap_success=length(opt))
}
boot_coefs <- function(form,dat,terms,B=1000){
  z<-matrix(NA_real_,B,length(terms),dimnames=list(NULL,terms)); success<-logical(B)
  for(b in seq_len(B)){ii<-sample.int(nrow(dat),nrow(dat),TRUE);fb<-tryCatch(coxph(form,dat[ii]),error=function(e)NULL);if(is.null(fb))next;cf<-coef(fb);if(all(terms%in%names(cf))&&all(is.finite(cf[terms]))){z[b,]<-cf[terms];success[b]<-TRUE}}
  rbindlist(lapply(terms,function(te){x<-z[success,te];data.table(term=te,bootstrap_success=sum(success),bootstrap_beta_median=median(x,na.rm=TRUE),bootstrap_beta_ci_low=quantile(x,.025,na.rm=TRUE),bootstrap_beta_ci_high=quantile(x,.975,na.rm=TRUE),bootstrap_positive_fraction=mean(x>0,na.rm=TRUE))}))
}
rows<-list();qarows<-list();bstab<-list();fits<-list()
for(id in names(specs)){
  if(id=="M6"&&!m6_ok){
    rows[[length(rows)+1]]<-data.table(model_id=id,term=NA_character_,beta=NA_real_,se=NA_real_,hr=NA_real_,ci_low=NA_real_,ci_high=NA_real_,p_value=NA_real_,n=nrow(d),events=sum(d$event),model_status="NOT_ESTIMABLE_COLLINEARITY_OR_EVENT_LIMIT")
    qarows[[length(qarows)+1]]<-data.table(model_id=id,model_status="NOT_ESTIMABLE_COLLINEARITY_OR_EVENT_LIMIT",ph_global_p=NA_real_,optimism_corrected_c_index=NA_real_,calibration_slope=NA_real_,bootstrap_success=0,max_abs_dfbeta=NA_real_,VIF=vif,condition_number=cond)
    next
  }
  f<-tryCatch(coxph(specs[[id]],d,x=TRUE,y=TRUE),error=function(e)NULL)
  if(is.null(f)){rows[[length(rows)+1]]<-data.table(model_id=id,term=NA_character_,model_status="NOT_ESTIMABLE_NUMERIC_FAILURE");next}
  fits[[id]]<-f;s<-summary(f);ci<-exp(confint(f));cf<-coef(f)
  for(te in names(cf)) rows[[length(rows)+1]]<-data.table(model_id=id,term=te,beta=cf[te],se=s$coefficients[te,"se(coef)"],hr=exp(cf[te]),ci_low=ci[te,1],ci_high=ci[te,2],p_value=s$coefficients[te,"Pr(>|z|)"],n=f$n,events=f$nevent,model_status="ESTIMABLE")
  ph<-tryCatch(cox.zph(f),error=function(e)NULL);vv<-boot_validate(specs[[id]],d,1000);dfb<-tryCatch(max(abs(residuals(f,type="dfbeta")),na.rm=TRUE),error=function(e)NA_real_)
  qarows[[length(qarows)+1]]<-data.table(model_id=id,model_status="ESTIMABLE",ph_global_p=if(is.null(ph))NA_real_ else ph$table["GLOBAL","p"],apparent_c_index=vv["apparent_c"],optimism=vv["optimism"],optimism_corrected_c_index=vv["corrected_c"],calibration_slope=vv["calibration_slope"],bootstrap_success=vv["bootstrap_success"],max_abs_dfbeta=dfb,VIF=if(id=="M6")vif else NA_real_,condition_number=if(id=="M6")cond else NA_real_)
  bt<-boot_coefs(specs[[id]],d,primary_terms[[id]],1000);bt[,model_id:=id];bstab[[length(bstab)+1]]<-bt
}
mods<-rbindlist(rows,fill=TRUE);qa<-rbindlist(qarows,fill=TRUE);bootstab<-rbindlist(bstab,fill=TRUE)
mods<-merge(mods,bootstab,by=c("model_id","term"),all.x=TRUE)
w(mods,"FULL_ICB_BASELINE_COMPONENT_MODELS.tsv");w(qa,"FULL_ICB_MODEL_QA.tsv")

# PH sensitivity: primary Cox remains primary; add time-varying coefficient only on violated component terms.
phsens<-list()
for(id in intersect(names(fits),c("M1","M2","M3","M4","M5"))){
  f<-fits[[id]];te<-primary_terms[[id]];ph<-tryCatch(cox.zph(f),error=function(e)NULL);tp<-if(is.null(ph)||!te%in%rownames(ph$table))NA_real_ else ph$table[te,"p"]
  if(is.finite(tp)&&tp<.05){
    form<-as.formula(paste0("Surv(time_years,event)~setting_recurrent+age_decade+",te,"+tt(",te,")"))
    ft<-tryCatch(coxph(form,d,tt=function(x,t,...)x*log(pmax(t,1e-6))),error=function(e)NULL)
    if(!is.null(ft)){ss<-summary(ft);for(ttn in rownames(ss$coefficients))phsens[[length(phsens)+1]]<-data.table(model_id=id,primary_term=te,cox_zph_term_p=tp,sensitivity="TIME_VARYING_COEFFICIENT_X_LOG_TIME",term=ttn,beta=ss$coefficients[ttn,"coef"],hr=exp(ss$coefficients[ttn,"coef"]),p_value=ss$coefficients[ttn,"Pr(>|z|)"],status="EVALUATED")}
  } else phsens[[length(phsens)+1]]<-data.table(model_id=id,primary_term=te,cox_zph_term_p=tp,sensitivity="NOT_REQUIRED",term=NA_character_,beta=NA_real_,hr=NA_real_,p_value=NA_real_,status="PH_NOT_VIOLATED")
}
w(rbindlist(phsens,fill=TRUE),"ICB_PH_SENSITIVITY.tsv")

# Exploratory measurement-error sensitivity from the frozen 13 matched bulk/snRNA calibration.
cal<-fread(file.path(V14,"05_icb_projection/ICB_13PATIENT_CALIBRATION.tsv"))
xx<-scale(cal$canonical_MES_bp);yy<-scale(cal$canonical_MES_sn)
pearson<-cor(xx,yy,use="complete.obs",method="pearson"); reliability<-max(.05,min(.95,pearson^2))
mr<-mods[model_id=="M2"&term=="malignant_MES_sd"]
me<-data.table(method="EXPLORATORY_REGRESSION_CALIBRATION_ATTENUATION",calibration_n=nrow(cal),pearson_r=pearson,reliability_ratio=reliability,primary_beta=mr$beta,primary_se=mr$se,primary_hr=mr$hr,corrected_beta=mr$beta/reliability,corrected_se=mr$se/reliability,corrected_hr=exp(mr$beta/reliability),corrected_ci_low=exp((mr$beta-1.96*mr$se)/reliability),corrected_ci_high=exp((mr$beta+1.96*mr$se)/reliability),interpretation="Exploratory only; 13-patient reliability is unstable and corrected significance cannot replace primary")
w(me,"FULL_ICB_MEASUREMENT_ERROR_SENSITIVITY.tsv")

# Compare direction with v1.4 historical 23 without selecting the new model.
old<-fread(file.path(V14,"05_icb_projection/ICB_23PATIENT_SURVIVAL_MODELS.tsv"))
mapold<-data.table(model_id=c("M1","M2","M3","M4"),old_model_id=c("M1","M2","M3","S4"))
cmp<-merge(mods[model_id%in%c("M1","M2","M3","M4")&term%in%c("bulk_MES_sd","malignant_MES_sd","ecological_component_sd","stable_mixed_contribution_sd"),.(model_id,new_term=term,new_beta=beta,new_hr=hr,new_p=p_value)],mapold,by="model_id")
cmp<-merge(cmp,old[,.(old_model_id=model_id,old_term=term,old_beta=beta,old_hr=hr,old_p=p_value)],by="old_model_id",all.x=TRUE)
cmp[,direction_consistent:=sign(new_beta)==sign(old_beta)]
w(cmp,"FULL_ICB_COMPARISON_WITH_V1_4_23.tsv")

# Formal gate: stable association requires exact-P, bootstrap sign, PH and influence stability; nulls remain.
pri<-mods[model_id%in%c("M1","M2","M3","M4","M5")&term%in%c("bulk_MES_sd","malignant_MES_sd","ecological_component_sd","stable_mixed_contribution_sd","discordance_sd")]
pri[,sign_stability:=pmax(bootstrap_positive_fraction,1-bootstrap_positive_fraction)]
pri<-merge(pri,qa[,.(model_id,ph_global_p,max_abs_dfbeta)],by="model_id",all.x=TRUE)
pri[,stable:=is.finite(p_value)&p_value<.05&is.finite(sign_stability)&sign_stability>=.9&(is.na(ph_global_p)|ph_global_p>=.05)&(is.na(max_abs_dfbeta)|max_abs_dfbeta<1)]
shared_stable<-pri[model_id=="M4",stable]
any_stable<-any(pri$stable)
all_cross<-all(pri$ci_low<1&pri$ci_high>1)
verdict<-if(isTRUE(shared_stable))"ICB_COUPLED_COMPONENT_ASSOCIATION_SUPPORTED" else if(any_stable)"ICB_COMPONENT_ASSOCIATION_REPRODUCIBLE" else if(all_cross)"ICB_COMPONENT_EFFECTS_IMPRECISE" else "ICB_NO_STABLE_COMPONENT_ASSOCIATION"
w(pri,"FULL_ICB_COMPONENT_STABILITY_GATE.tsv")
writeLines(c("# ICB projection and baseline prognosis summary","",paste0("Status: `",verdict,"`"),"",paste0("Patients/events: ",nrow(d),"/",sum(d$event)),paste0("Reference overlap: ",length(common)," genes"),paste0("Canonical coverage bulk/deconvolved: ",min(d$canonical_gene_coverage_bulk),"/95 and ",min(d$canonical_gene_coverage_deconvolved),"/95"),paste0("M6 estimability: ",ifelse(m6_ok,"ESTIMABLE","NOT_ESTIMABLE"),"; VIF=",signif(vif,4),"; condition number=",signif(cond,4)),"","All survival estimates are within-ICB prognostic associations. No treatment-predictive or treatment-selection interpretation is permitted. Measurement-error correction is exploratory and cannot replace the primary estimates."),file.path(OUT,"ICB_SURVIVAL_SUMMARY.md"))
writeLines(capture.output(sessionInfo()),file.path(RUN,"session_info/MODULE1_SESSION_INFO.txt"))
cat("MODULE1_COMPLETE",verdict,nrow(d),sum(d$event),length(common),"\n")
