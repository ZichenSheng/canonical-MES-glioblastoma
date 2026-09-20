#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(survival)})
set.seed(16061)
DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset="data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset="results")
REPO_ROOT <- Sys.getenv("GBM_MES_REPO_ROOT", unset=".")
RUN<-file.path(RESULT_ROOT,"clinical_glass")
V15<-file.path(RESULT_ROOT,"clinical_icb")
OUT<-file.path(RUN,"06_icb_drug_sensitivity");dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
d<-fread(file.path(V15,"FINAL_FULL_ICB_COHORT_V1_5.tsv"),na.strings=c("NA",""))
stopifnot(nrow(d)==40,uniqueN(d$participant_id)==40)
d[,agent_source:=fifelse(!is.na(`Checkpoint blockade immuno`)&nzchar(`Checkpoint blockade immuno`),`Checkpoint blockade immuno`,icb_drug)]
d[,agent_norm:=tolower(trimws(agent_source))]
d[,agent_standard:=fcase(grepl("multiple|\\+|,",agent_norm),"combination",grepl("nivolumab",agent_norm),"nivolumab",grepl("pembrolizumab",agent_norm),"pembrolizumab",grepl("durvalumab",agent_norm),"durvalumab",default="unknown")]
d[,target_class:=fcase(agent_standard%in%c("nivolumab","pembrolizumab"),"anti-PD-1",agent_standard=="durvalumab","anti-PD-L1",agent_standard=="combination","combination/other",default="unknown")]
d[,anti_pd1_only:=agent_standard%in%c("nivolumab","pembrolizumab")]
dist<-d[,.(n=.N,events=sum(event),newly_diagnosed=sum(setting_recurrent==0),recurrent=sum(setting_recurrent==1),dexamethasone_missing=sum(is.na(`Dex at start of ICB_score`)),bevacizumab_during_missing=sum(is.na(`Bevacizumab use during ICB_score`)),concomitant_drug_summary=paste(sort(unique(na.omit(as.character(Chemotherapy)))),collapse=";"),source_cohort=paste(sort(unique(na.omit(as.character(`Dx at time of ICB (DFCI)`)))),collapse=";")),by=.(agent_standard,target_class)]
fwrite(dist,file.path(OUT,"ICB_AGENT_DISTRIBUTION.tsv"),sep="\t",quote=FALSE,na="NA")

fit_score<-function(dd,score,analysis_set,B=1000L){
 f<-as.formula(paste0("Surv(time_years,event)~setting_recurrent+age_decade+",score));fit<-coxph(f,data=dd,ties="efron",x=TRUE,y=TRUE)
 sm<-summary(fit);j<-match(score,rownames(sm$coefficients));z<-cox.zph(fit,transform="km")
 boot<-rep(NA_real_,B);for(b in seq_len(B)){ix<-sample.int(nrow(dd),nrow(dd),TRUE);boot[b]<-tryCatch(coef(coxph(f,data=dd[ix],ties="efron"))[score],error=function(e)NA_real_)};boot<-boot[is.finite(boot)]
 db<-residuals(fit,type="dfbeta");if(is.null(dim(db)))db<-matrix(db,ncol=1,dimnames=list(NULL,names(coef(fit))));if(is.null(colnames(db)))colnames(db)<-names(coef(fit));jj<-match(score,colnames(db))
  dbv<-if(is.na(jj))numeric()else db[,jj];maxdb<-if(any(is.finite(dbv)))max(abs(dbv[is.finite(dbv)]))else NA_real_
  data.table(analysis_set=analysis_set,score=score,beta=sm$coefficients[j,"coef"],se=sm$coefficients[j,"se(coef)"],hr=sm$conf.int[j,"exp(coef)"],ci_low=sm$conf.int[j,"lower .95"],ci_high=sm$conf.int[j,"upper .95"],p_value=sm$coefficients[j,"Pr(>|z|)"],n=nrow(dd),events=sum(dd$event),newly_diagnosed=sum(dd$setting_recurrent==0),recurrent=sum(dd$setting_recurrent==1),ph_score_p=z$table[score,"p"],ph_global_p=z$table["GLOBAL","p"],bootstrap_success=length(boot),bootstrap_beta_median=median(boot),bootstrap_ci_low=quantile(boot,.025),bootstrap_ci_high=quantile(boot,.975),bootstrap_negative_fraction=mean(boot<0),bootstrap_sign_agreement=mean(sign(boot)==sign(sm$coefficients[j,"coef"])),max_abs_dfbeta=maxdb)
}
scores<-c("bulk_MES_sd","stable_mixed_contribution_sd")
res<-rbindlist(lapply(scores,function(s)rbind(fit_score(d,s,"FULL_ICB_40"),fit_score(d[anti_pd1_only==TRUE],s,"ANTI_PD1_ONLY"))))
full<-res[analysis_set=="FULL_ICB_40",.(score,full_beta=beta,full_ci_low=log(ci_low),full_ci_high=log(ci_high))]
res<-merge(res,full,by="score",all.x=TRUE)
res[,beta_shift_from_full:=beta-full_beta]
res[,direction_status:=fcase(analysis_set=="FULL_ICB_40","REFERENCE",beta<0&beta>=full_ci_low&beta<=full_ci_high&abs(beta_shift_from_full)<=.35,"STABLE",beta<0,"DIRECTIONALLY_STABLE",n>=15&events>=10,"AGENT_SENSITIVE",default="UNDERPOWERED")]
res[,ph_status:=fifelse(is.finite(ph_score_p)&ph_score_p>=.05&is.finite(ph_global_p)&ph_global_p>=.05,"PASS","PH_SENSITIVE")]
fwrite(res,file.path(OUT,"ICB_ANTI_PD1_SENSITIVITY.tsv"),sep="\t",quote=FALSE,na="NA")

loo<-list();q<-0
for(agent in c("nivolumab","pembrolizumab","durvalumab")){
 ex<-d[agent_standard==agent];keep<-d[agent_standard!=agent];gate=nrow(ex)>=8&&sum(ex$event)>=5&&nrow(keep)>=8&&sum(keep$event)>=5
 for(s in scores){q<-q+1;if(gate){r<-fit_score(keep,s,paste0("LEAVE_",toupper(agent),"_OUT"),B=1000);r[,`:=`(excluded_agent=agent,excluded_n=nrow(ex),excluded_events=sum(ex$event),execution_gate="PASS_EXECUTED")];loo[[q]]<-r}else loo[[q]]<-data.table(analysis_set=paste0("LEAVE_",toupper(agent),"_OUT"),score=s,beta=NA_real_,se=NA_real_,hr=NA_real_,ci_low=NA_real_,ci_high=NA_real_,p_value=NA_real_,n=nrow(keep),events=sum(keep$event),newly_diagnosed=sum(keep$setting_recurrent==0),recurrent=sum(keep$setting_recurrent==1),ph_score_p=NA_real_,ph_global_p=NA_real_,bootstrap_success=0L,bootstrap_beta_median=NA_real_,bootstrap_ci_low=NA_real_,bootstrap_ci_high=NA_real_,bootstrap_negative_fraction=NA_real_,bootstrap_sign_agreement=NA_real_,max_abs_dfbeta=NA_real_,excluded_agent=agent,excluded_n=nrow(ex),excluded_events=sum(ex$event),execution_gate="NOT_EXECUTED_MIN_8_PATIENTS_5_EVENTS")
 }
}
loo<-rbindlist(loo,fill=TRUE);loo<-merge(loo,full,by="score",all.x=TRUE);loo[,beta_shift_from_full:=beta-full_beta]
fwrite(loo,file.path(OUT,"ICB_LEAVE_ONE_AGENT_OUT.tsv"),sep="\t",quote=FALSE,na="NA")
baseline<-fread(file.path(V15,"FINAL_ICB_BASELINE_COMPONENT_MODELS_V1_5.tsv"));checks<-merge(res[analysis_set=="FULL_ICB_40",.(score,beta,hr)],baseline[term%in%scores,.(score=term,v15_beta=beta,v15_hr=hr)],by="score");checks[,pass:=abs(beta-v15_beta)<1e-12&abs(hr-v15_hr)<1e-12]
if(!all(checks$pass))stop("full ICB model did not reproduce exactly")
ap<-res[analysis_set=="ANTI_PD1_ONLY"]
gate<-if(all(ap$direction_status=="STABLE")&&all(ap$ph_status=="PASS"))"ICB_SHARED_ASSOCIATION_STABLE_IN_ANTI_PD1" else if(all(ap$beta<0))"ICB_SHARED_ASSOCIATION_DIRECTIONALLY_STABLE_IN_ANTI_PD1" else if(min(ap$n)<15||min(ap$events)<10)"ICB_DRUG_SENSITIVITY_UNDERPOWERED" else "ICB_SHARED_ASSOCIATION_AGENT_SENSITIVE"
writeLines(c("# ICB drug sensitivity gate","",paste0("Status: `",gate,"`"),"",paste0("Frozen cohort: 40 patients / ",sum(d$event)," events"),paste0("Anti-PD-1-only: ",sum(d$anti_pd1_only)," patients / ",sum(d[anti_pd1_only==TRUE]$event)," events"),"Full-cohort coefficients reproduce v1.5 exactly. Both effects retain direction and similar magnitude in anti-PD-1-only, but the anti-PD-1 clinical model has a global PH violation; the formal claim is therefore directional rather than fully stable. Anti-PD-1 significance is not a gate. Leave-durvalumab-out is not executed because the excluded group has fewer than eight patients."),file.path(OUT,"ICB_DRUG_SENSITIVITY_GATE.md"))
fwrite(checks,file.path(OUT,"ICB_FULL_MODEL_REPRODUCTION_QA.tsv"),sep="\t",quote=FALSE,na="NA")
writeLines(capture.output(sessionInfo()),file.path(RUN,"session_info/ICB_DRUG_SENSITIVITY_SESSION_INFO.txt"))
cat(paste0('{"status":"',gate,'","anti_pd1_n":',sum(d$anti_pd1_only),',"anti_pd1_events":',sum(d[anti_pd1_only==TRUE]$event),'}\n'))

# Randomized post-one-dose versus pre-exposure tissue contrast (GSE121810).
PD1_OUT <- file.path(RUN,"07_randomized_pd1_perturbation");dir.create(PD1_OUT,recursive=TRUE,showWarnings=FALSE)
pd1 <- fread(file.path(DATA_ROOT,"clinical","gse121810","NEOADJUVANT_PD1_SAMPLE_SCORES.tsv"))[study=="GSE121810"]
stopifnot(nrow(pd1)==29L,pd1[arm=="neoadjuvant",.N]==14L,pd1[arm=="adjuvant_only",.N]==15L,
          uniqueN(pd1$patient_id)==29L,all(pd1$patient_independent))
pd1[,neo:=as.integer(arm=="neoadjuvant")]

hedges_g <- function(y,g){
 y1<-y[g==1L];y0<-y[g==0L];n1<-length(y1);n0<-length(y0)
 sp<-sqrt(((n1-1)*var(y1)+(n0-1)*var(y0))/(n1+n0-2));d0<-(mean(y1)-mean(y0))/sp
 (1-3/(4*(n1+n0)-9))*d0
}
randomized_arm_test <- function(var,family,B=2000L,P=10000L,seed=1L){
 y<-pd1[[var]];ok<-is.finite(y);y<-y[ok];g<-pd1$neo[ok]
 obs_mean<-mean(y[g==1L])-mean(y[g==0L]);obs_median<-median(y[g==1L])-median(y[g==0L]);obs_g<-hedges_g(y,g)
 fit<-lm(y~g);set.seed(seed)
 perm<-vapply(seq_len(P),function(i){gp<-sample(g,replace=FALSE);mean(y[gp==1L])-mean(y[gp==0L])},numeric(1))
 p_perm<-(1+sum(abs(perm)>=abs(obs_mean)))/(P+1);set.seed(seed+1L)
 boot<-replicate(B,{
  i1<-sample(which(g==1L),sum(g==1L),replace=TRUE);i0<-sample(which(g==0L),sum(g==0L),replace=TRUE)
  ii<-c(i1,i0);gb<-c(rep(1L,length(i1)),rep(0L,length(i0)));yb<-y[ii]
  c(mean_diff=mean(yb[gb==1L])-mean(yb[gb==0L]),median_diff=median(yb[gb==1L])-median(yb[gb==0L]),hedges_g=hedges_g(yb,gb))
 })
 data.table(endpoint=var,family=family,n=length(y),neoadjuvant_n=sum(g==1L),adjuvant_only_n=sum(g==0L),
  adjuvant_mean=mean(y[g==0L]),neoadjuvant_mean=mean(y[g==1L]),mean_difference_neoadjuvant_minus_adjuvant=obs_mean,
  mean_difference_ci_low=quantile(boot["mean_diff",],.025,na.rm=TRUE),mean_difference_ci_high=quantile(boot["mean_diff",],.975,na.rm=TRUE),
  adjuvant_median=median(y[g==0L]),neoadjuvant_median=median(y[g==1L]),median_difference_neoadjuvant_minus_adjuvant=obs_median,
  median_difference_ci_low=quantile(boot["median_diff",],.025,na.rm=TRUE),median_difference_ci_high=quantile(boot["median_diff",],.975,na.rm=TRUE),
  standardized_mean_difference_hedges_g=obs_g,smd_ci_low=quantile(boot["hedges_g",],.025,na.rm=TRUE),
  smd_ci_high=quantile(boot["hedges_g",],.975,na.rm=TRUE),linear_model_p=summary(fit)$coefficients["g","Pr(>|t|)"],
  permutation_p=p_perm,permutation_n=P,bootstrap_n=B)
}
components<-c("bulk_MES","stable_mixed_contribution","malignant_biased_contribution","ecology_biased_contribution",
              "malignant_MES","ecological_component","discordance")
immune<-c("IFN_response","T_cell_fraction","myeloid_fraction","HLA_I","cell_cycle")
pd1_rows<-rbindlist(lapply(seq_along(components),function(i)randomized_arm_test(components[i],"PRESPECIFIED_MES_COMPONENT_FAMILY",seed=2026081206L+i*100L)))
pd1_rows[,BH_FDR:=p.adjust(permutation_p,method="BH")]
pd1_rows[,status:="EVALUATED_RANDOMIZED_PHARMACODYNAMIC_CONTRAST"]
pd1_rows[,interpretation_boundary:="post-one-dose versus pre-exposure randomized resection-tissue contrast; not baseline prediction"]
fwrite(pd1_rows,file.path(PD1_OUT,"GSE121810_RANDOMIZED_PD1_PERTURBATION.tsv"),sep="\t",quote=FALSE,na="NA")
immune_rows<-rbindlist(lapply(seq_along(immune),function(i)randomized_arm_test(immune[i],"SUPPORTIVE_IMMUNE_PROGRAM",seed=2026081606L+i*100L)))
immune_rows[,BH_FDR_supportive:=p.adjust(permutation_p,method="BH")]
fwrite(immune_rows,file.path(PD1_OUT,"GSE121810_IMMUNE_PROGRAM_ARM_CONTRAST.tsv"),sep="\t",quote=FALSE,na="NA")
pd1_gate<-if(any(pd1_rows$BH_FDR<.10))"RANDOMIZED_PHARMACODYNAMIC_RESPONSIVENESS_SUPPORTED" else "NO_REPRODUCIBLE_PD1_PERTURBATION_OF_MES_ARCHITECTURE"
writeLines(pd1_gate,file.path(PD1_OUT,"GSE121810_PERTURBATION_STATUS.txt"))
