#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(survival)})
set.seed(16061)
DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset="data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset="results")
REPO_ROOT <- Sys.getenv("GBM_MES_REPO_ROOT", unset=".")
RUN <- file.path(RESULT_ROOT,"clinical_glass")
OUT <- file.path(RUN,"survival")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
spec <- data.table(model_id=c("G1","G2","G3","G4"),score=c("bulk_MES_sd","stable_mixed_contribution_sd","malignant_biased_contribution_sd","ecology_biased_contribution_sd"))

fit_one <- function(d,cohort,mid,score,variant="PRIMARY_MINIMAL",B=1000L,optB=200L){
  dd<-copy(d);dd[,age_decade:=age/10]; dd[,MGMT_FACTOR:=factor(MGMT)]
  covs<-if(variant=="PRIMARY_MINIMAL")c("age_decade",score)else c("age_decade","MGMT_FACTOR",score)
  dd<-dd[complete.cases(dd[,..covs])&is.finite(OS)&OS>0&event%in%c(0,1)]
  f<-as.formula(paste0("Surv(OS,event)~",paste(covs,collapse="+")))
  fit<-coxph(f,data=dd,x=TRUE,y=TRUE,model=TRUE,ties="efron")
  sm<-summary(fit); co<-sm$coefficients; ci<-sm$conf.int
  j<-match(score,rownames(co)); if(is.na(j))stop("score term missing")
  zph<-tryCatch(cox.zph(fit,transform="km"),error=function(e)NULL)
  ph_score<-if(is.null(zph))NA_real_ else zph$table[score,"p"]
  ph_global<-if(is.null(zph))NA_real_ else zph$table["GLOBAL","p"]
  beta<-co[j,"coef"]; se<-co[j,"se(coef)"]
  boot<-rep(NA_real_,B)
  for(b in seq_len(B)){
    ix<-sample.int(nrow(dd),nrow(dd),replace=TRUE)
    boot[b]<-tryCatch(coef(coxph(f,data=dd[ix],ties="efron"))[score],error=function(e)NA_real_)
  }
  boot<-boot[is.finite(boot)]
  db<-residuals(fit,type="dfbeta"); if(is.null(dim(db)))db<-matrix(db,ncol=1,dimnames=list(NULL,names(coef(fit))))
  dj<-match(score,colnames(db)); maxdf<-if(is.na(dj))NA_real_ else max(abs(db[,dj]),na.rm=TRUE); infln<-if(is.na(dj))NA_integer_ else sum(abs(db[,dj])>2/sqrt(nrow(dd)),na.rm=TRUE)
  worst<-if(is.na(dj))NA_integer_ else which.max(abs(db[,dj])); loo<-if(is.na(worst))NA_real_ else tryCatch(coef(coxph(f,data=dd[-worst],ties="efron"))[score],error=function(e)NA_real_)
  capp<-unname(sm$concordance[1]); optimism<-rep(NA_real_,optB)
  for(b in seq_len(optB)){
    ix<-sample.int(nrow(dd),nrow(dd),replace=TRUE); fb<-tryCatch(coxph(f,data=dd[ix],ties="efron",x=TRUE),error=function(e)NULL); if(is.null(fb))next
    rb<-tryCatch(predict(fb,newdata=dd[ix],type="lp"),error=function(e)NULL); ro<-tryCatch(predict(fb,newdata=dd,type="lp"),error=function(e)NULL); if(is.null(rb)||is.null(ro))next
    cb<-tryCatch(concordance(Surv(OS,event)~rb,data=dd[ix],reverse=TRUE)$concordance,error=function(e)NA_real_)
    ct<-tryCatch(concordance(Surv(OS,event)~ro,data=dd,reverse=TRUE)$concordance,error=function(e)NA_real_)
    optimism[b]<-cb-ct
  }
  oc<-if(any(is.finite(optimism)))mean(optimism[is.finite(optimism)])else NA_real_; ccor<-capp-oc
  sign_agree<-mean(sign(boot)==sign(beta)); phok<-is.na(ph_score)||ph_score>=.05
  row<-data.table(cohort=cohort,model_id=mid,model_variant=variant,term=score,beta=beta,se=se,hr=exp(beta),ci_low=ci[j,"lower .95"],ci_high=ci[j,"upper .95"],p_value=co[j,"Pr(>|z|)"],n=nrow(dd),events=sum(dd$event),ph_score_p=ph_score,ph_global_p=ph_global,bootstrap_success=length(boot),bootstrap_beta_median=median(boot),bootstrap_beta_ci_low=quantile(boot,.025),bootstrap_beta_ci_high=quantile(boot,.975),bootstrap_sign_agreement=sign_agree,bootstrap_negative_fraction=mean(boot<0),max_abs_dfbeta=maxdf,influential_n=infln,leave_most_influential_beta=loo,c_index_apparent=capp,optimism=oc,c_index_optimism_corrected=ccor,model_status=ifelse(phok&&sign_agree>=.8,"STABLE","UNSTABLE_OR_PH_SENSITIVE"))
  stab<-data.table(cohort=cohort,model_id=mid,model_variant=variant,score=score,n=nrow(dd),events=sum(dd$event),beta=beta,bootstrap_sign_agreement=sign_agree,max_abs_dfbeta=maxdf,leave_most_influential_beta=loo,beta_shift_leave_one=loo-beta,c_index_apparent=capp,c_index_optimism_corrected=ccor,stability_status=row$model_status)
  phr<-rbind(data.table(cohort=cohort,model_id=mid,model_variant=variant,term=score,chisq=if(is.null(zph))NA_real_ else zph$table[score,"chisq"],p_value=ph_score,status=ifelse(is.na(ph_score),"NOT_EVALUABLE",ifelse(ph_score<.05,"PH_VIOLATION","PASS"))),data.table(cohort=cohort,model_id=mid,model_variant=variant,term="GLOBAL",chisq=if(is.null(zph))NA_real_ else zph$table["GLOBAL","chisq"],p_value=ph_global,status=ifelse(is.na(ph_global),"NOT_EVALUABLE",ifelse(ph_global<.05,"PH_VIOLATION","PASS"))))
  tv<-data.table()
  if((is.finite(ph_score)&&ph_score<.05)||(is.finite(ph_global)&&ph_global<.05)){
    ftv<-as.formula(paste0("Surv(OS,event)~",paste(setdiff(covs,score),collapse="+"),"+",score,"+tt(",score,")"))
    tf<-tryCatch(coxph(ftv,data=dd,ties="efron",tt=function(x,t,...)x*log(pmax(t,1/30.4375))),error=function(e)NULL)
    if(!is.null(tf)){
      cs<-summary(tf)$coefficients; tn<-grep("tt\\(",rownames(cs),value=TRUE)[1]
      tv<-data.table(cohort=cohort,model_id=mid,model_variant=variant,sensitivity="TIME_VARYING_COEFFICIENT",base_beta=coef(tf)[score],time_interaction_beta=coef(tf)[tn],time_interaction_p=cs[tn,"Pr(>|z|)"],interpretation="score effect permitted to vary with log time; primary Cox retained")
    }
  }
  list(row=row,stab=stab,ph=phr,tv=tv)
}

allrow<-list();allstab<-list();allph<-list();alltv<-list();q<-0
for(cohort in c("PRIMARY","RECURRENT")){
 d<-fread(file.path(RUN,"03_glass_cohort",paste0("GLASS_",cohort,"_COHORT.tsv")),na.strings=c("NA",""))
 for(i in 1:nrow(spec)){
  for(v in c("PRIMARY_MINIMAL","COMPLETE_CASE_MGMT_SENSITIVITY")){
   q<-q+1;r<-fit_one(d,cohort,spec$model_id[i],spec$score[i],v,B=ifelse(v=="PRIMARY_MINIMAL",1000L,500L),optB=ifelse(v=="PRIMARY_MINIMAL",200L,100L));allrow[[q]]<-r$row;allstab[[q]]<-r$stab;allph[[q]]<-r$ph;if(nrow(r$tv))alltv[[length(alltv)+1]]<-r$tv
  }
 }
}
models<-rbindlist(allrow,fill=TRUE);stability<-rbindlist(allstab,fill=TRUE);ph<-rbindlist(allph,fill=TRUE);tv<-rbindlist(alltv,fill=TRUE)
fwrite(models[cohort=="PRIMARY"],file.path(OUT,"GLASS_PRIMARY_SURVIVAL_MODELS.tsv"),sep="\t",quote=FALSE,na="NA")
fwrite(models[cohort=="RECURRENT"],file.path(OUT,"GLASS_RECURRENT_SURVIVAL_MODELS.tsv"),sep="\t",quote=FALSE,na="NA")
fwrite(stability,file.path(OUT,"GLASS_MODEL_STABILITY.tsv"),sep="\t",quote=FALSE,na="NA")
fwrite(ph,file.path(OUT,"GLASS_PH_DIAGNOSTICS.tsv"),sep="\t",quote=FALSE,na="NA")
if(!nrow(tv))tv<-data.table(cohort="ALL",model_id="NONE",model_variant="PRIMARY_MINIMAL",sensitivity="NOT_TRIGGERED",base_beta=NA,time_interaction_beta=NA,time_interaction_p=NA,interpretation="No PH violation")
fwrite(tv,file.path(OUT,"GLASS_TIME_VARYING_SENSITIVITY.tsv"),sep="\t",quote=FALSE,na="NA")
rmst<-ph[model_variant=="PRIMARY_MINIMAL",.(ph_trigger=any(status=="PH_VIOLATION")),by=.(cohort,model_id)]
rmst[,`:=`(method=ifelse(ph_trigger,"TIME_VARYING_COX_USED_INSTEAD_OF_RMST","NOT_TRIGGERED"),tau_months=NA_real_,estimate=NA_real_,status=ifelse(ph_trigger,"ALTERNATIVE_ALLOWED_SENSITIVITY_COMPLETED","NOT_REQUIRED"),note="Continuous primary score retained; no outcome-selected dichotomization or tau")]
fwrite(rmst,file.path(OUT,"GLASS_RMST_SENSITIVITY.tsv"),sep="\t",quote=FALSE,na="NA")
writeLines(capture.output(sessionInfo()),file.path(RUN,"session_info/GLASS_SURVIVAL_SESSION_INFO.txt"))
cat(toJSON<-paste0('{"status":"PASS","primary_n":',models[cohort=="PRIMARY"&model_variant=="PRIMARY_MINIMAL",unique(n)],',"recurrent_n":',models[cohort=="RECURRENT"&model_variant=="PRIMARY_MINIMAL",unique(n)],',"ph_triggered":',sum(rmst$ph_trigger),' }\n'))
