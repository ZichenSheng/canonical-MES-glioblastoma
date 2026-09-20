#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(survival); library(splines)})

DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset="data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset="results")
REPO_ROOT <- Sys.getenv("GBM_MES_REPO_ROOT", unset=".")
run <- file.path(RESULT_ROOT,"clinical_context")
seed <- 17032026L
B <- 500L
set.seed(seed)

readt <- function(p) read.delim(file.path(run,p), check.names=FALSE, stringsAsFactors=FALSE)
writet <- function(x,p) write.table(x,file.path(run,p),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
icb <- readt("03_score_scale_harmonization/ICB_HARMONIZED_SCORES.tsv")
glass <- readt("03_score_scale_harmonization/GLASS_HARMONIZED_SCORES.tsv")
dall <- rbind(icb,glass)
dall$age_z <- as.numeric(scale(dall$age))

split_month <- function(d) {
  out <- vector("list",nrow(d))
  for(i in seq_len(nrow(d))) {
    tt <- max(as.numeric(d$time_months[i]),0.05)
    cuts <- unique(c(seq(0,floor(tt),by=1),tt)); cuts <- sort(cuts[cuts<=tt])
    if(tail(cuts,1)<tt) cuts <- c(cuts,tt)
    if(length(cuts)<2) cuts <- c(0,tt)
    start <- head(cuts,-1); stop <- tail(cuts,-1); pt <- stop-start
    ev <- rep(0,length(pt)); if(as.numeric(d$event[i])==1) ev[length(ev)] <- 1
    out[[i]] <- data.frame(patient_id=d$patient_id[i],start=start,stop=stop,pt=pt,log_mid=log((start+stop)/2+0.05),event_interval=ev,score=d$score[i],age_z=d$age_z[i])
  }
  do.call(rbind,out)
}

fit_one <- function(d,score_col) {
  d$score <- d[[score_col]]
  # Re-standardize age within context; predictions are at mean age.
  d$age_z <- if(sd(d$age,na.rm=TRUE)>0) (d$age-mean(d$age,na.rm=TRUE))/sd(d$age,na.rm=TRUE) else 0
  sp <- split_month(d)
  fit <- try(glm(event_interval ~ ns(log_mid,df=3) + age_z + score + score:ns(log_mid,df=2),
                 offset=log(pt),family=poisson(),data=sp,control=glm.control(maxit=80)),silent=TRUE)
  if(inherits(fit,"try-error") || !isTRUE(fit$converged)) return(NULL)
  list(fit=fit,n=nrow(d),events=sum(d$event),rank=fit$rank,ncoef=length(coef(fit)))
}

predict_curve <- function(obj,score,tau,dt=0.05) {
  grid <- seq(dt/2,tau-dt/2,by=dt)
  nd <- data.frame(log_mid=log(grid+0.05),age_z=0,score=score,pt=dt)
  haz <- exp(as.numeric(predict(obj$fit,newdata=nd,type="link")) - log(dt))
  surv <- exp(-cumsum(haz*dt))
  list(time=grid,surv=surv,rmst=sum(surv*dt),haz=haz)
}

effects <- function(obj,tau,times=c(12,18)) {
  lo <- predict_curve(obj,0.25,max(c(tau,times)))
  hi <- predict_curve(obj,0.75,max(c(tau,times)))
  rmst_lo <- sum(lo$surv[lo$time<=tau]*0.05); rmst_hi <- sum(hi$surv[hi$time<=tau]*0.05)
  fixed <- sapply(times,function(tt){j<-which.min(abs(lo$time-tt)); hi$surv[j]-lo$surv[j]})
  list(rmst=rmst_hi-rmst_lo,fixed=fixed)
}

loghr_time <- function(obj,times) {
  out <- numeric(length(times))
  for(i in seq_along(times)) {
    tt <- times[i]
    nd <- data.frame(log_mid=rep(log(tt+0.05),2),age_z=0,score=c(.25,.75),pt=1)
    lp <- predict(obj$fit,newdata=nd,type="link")
    out[i] <- as.numeric(lp[2]-lp[1])
  }
  out
}

boot_pair <- function(di,dg,score_col,tau,times_grid) {
  store <- matrix(NA_real_,nrow=B,ncol=3+length(times_grid))
  colnames(store) <- c("rmst_delta","fixed12_delta","fixed18_delta",paste0("tv",times_grid))
  for(b in seq_len(B)) {
    bi <- di[sample.int(nrow(di),replace=TRUE),,drop=FALSE]
    bg <- dg[sample.int(nrow(dg),replace=TRUE),,drop=FALSE]
    # Entire score transformation is repeated within each bootstrap cohort.
    ri <- rank(bi[[sub("_percentile$","_raw",score_col)]],ties.method="average")
    rg <- rank(bg[[sub("_percentile$","_raw",score_col)]],ties.method="average")
    bi[[score_col]] <- (ri-.5)/nrow(bi); bg[[score_col]] <- (rg-.5)/nrow(bg)
    fi <- fit_one(bi,score_col); fg <- fit_one(bg,score_col)
    if(is.null(fi)||is.null(fg)) next
    ei <- effects(fi,tau); eg <- effects(fg,tau)
    store[b,1] <- ei$rmst-eg$rmst
    store[b,2:3] <- ei$fixed-eg$fixed
    store[b,4:ncol(store)] <- loghr_time(fi,times_grid)-loghr_time(fg,times_grid)
  }
  store
}

qci <- function(x) {x<-x[is.finite(x)]; if(length(x)<20) return(c(NA,NA,NA,NA)); c(quantile(x,.025,names=FALSE),quantile(x,.975,names=FALSE),2*min(mean(x<=0),mean(x>=0)),length(x))}

run_main <- function() {
rmst_rows <- list(); fixed_rows <- list(); tv_rows <- list(); coef_rows <- list(); qa_rows <- list(); cox_rows <- list(); k<-1
times_grid <- c(1,3,6,9,12,18,24)
for(stage in c("newly_diagnosed","recurrent")) {
  tau <- if(stage=="newly_diagnosed") 18 else 12
  di <- subset(dall,stage==stage & context=="ICB")
  dg <- subset(dall,stage==stage & context=="GLASS")
  # Avoid R's subset name collision.
  di <- dall[dall$stage==stage & dall$context=="ICB",,drop=FALSE]
  dg <- dall[dall$stage==stage & dall$context=="GLASS",,drop=FALSE]
  for(bio in c("stable_mixed","bulk_MES")) {
    sc <- paste0(bio,"_percentile")
    fi <- fit_one(di,sc); fg <- fit_one(dg,sc)
    if(is.null(fi)||is.null(fg)) next
    ei <- effects(fi,tau); eg <- effects(fg,tau)
    boot <- boot_pair(di,dg,sc,tau,times_grid)
    ci <- qci(boot[,"rmst_delta"])
    rmst_rows[[length(rmst_rows)+1]] <- data.frame(stage=stage,biomarker=bio,scale="percentile_IQR_0.75_vs_0.25",tau_months=tau,icb_effect_months=ei$rmst,glass_effect_months=eg$rmst,rmst_interaction_months=ei$rmst-eg$rmst,ci_low=ci[1],ci_high=ci[2],bootstrap_p=ci[3],bootstrap_success=ci[4],icb_n=nrow(di),icb_events=sum(di$event),glass_n=nrow(dg),glass_events=sum(dg$event),interpretation="observational origin-mismatched context contrast")
    for(j in 1:2) {
      tt <- c(12,18)[j]; cj <- qci(boot[,c("fixed12_delta","fixed18_delta")[j]])
      fixed_rows[[length(fixed_rows)+1]] <- data.frame(stage=stage,biomarker=bio,scale="percentile_IQR_0.75_vs_0.25",time_months=tt,icb_survival_probability_effect=ei$fixed[j],glass_survival_probability_effect=eg$fixed[j],context_risk_difference=ei$fixed[j]-eg$fixed[j],ci_low=cj[1],ci_high=cj[2],bootstrap_p=cj[3],bootstrap_success=cj[4],icb_at_risk=sum(di$time_months>=tt),glass_at_risk=sum(dg$time_months>=tt),support=ifelse(sum(di$time_months>=tt)>=5 & sum(dg$time_months>=tt)>=5,"ADEQUATE","LIMITED"))
    }
    lhi <- loghr_time(fi,times_grid); lhg <- loghr_time(fg,times_grid)
    for(j in seq_along(times_grid)) {
      cj <- qci(boot[,paste0("tv",times_grid[j])])
      tv_rows[[length(tv_rows)+1]] <- data.frame(stage=stage,biomarker=bio,time_months=times_grid[j],icb_logHR_IQR=lhi[j],glass_logHR_IQR=lhg[j],delta_logHR=lhi[j]-lhg[j],ci_low=cj[1],ci_high=cj[2],bootstrap_p=cj[3],bootstrap_success=cj[4],support=ifelse(sum(di$time_months>=times_grid[j])>=5 & sum(dg$time_months>=times_grid[j])>=5,"ADEQUATE","LIMITED"))
    }
    for(ctx in c("ICB","GLASS")) {
      obj <- if(ctx=="ICB") fi else fg; dd <- if(ctx=="ICB") di else dg
      cf <- coef(summary(obj$fit))
      for(term in rownames(cf)) coef_rows[[length(coef_rows)+1]] <- data.frame(stage=stage,context=ctx,biomarker=bio,model="piecewise_exponential_natural_spline",baseline_spline_df=3,time_varying_score_df=2,term=term,estimate=cf[term,1],se=cf[term,2],p_value=cf[term,4],n=nrow(dd),events=sum(dd$event),converged=obj$fit$converged,rank=obj$rank,ncoef=obj$ncoef)
      qa_rows[[length(qa_rows)+1]] <- data.frame(stage=stage,context=ctx,biomarker=bio,n=nrow(dd),events=sum(dd$event),tau_months=tau,at_risk_tau=sum(dd$time_months>=tau),baseline_spline_df=3,time_varying_score_df=2,interval_width_months=1,converged=obj$fit$converged,rank=obj$rank,ncoef=obj$ncoef,status=ifelse(obj$fit$converged & obj$rank==obj$ncoef,"PASS","LIMITED"))
      dd$age_z <- if(sd(dd$age,na.rm=TRUE)>0) (dd$age-mean(dd$age,na.rm=TRUE))/sd(dd$age,na.rm=TRUE) else 0
      fcox <- coxph(as.formula(paste0("Surv(time_months,event)~",sc,"+age_z")),data=dd,x=TRUE)
      cz <- cox.zph(fcox)
      rr <- summary(fcox)$coef[sc,]
      cox_rows[[length(cox_rows)+1]] <- data.frame(stage=stage,context=ctx,biomarker=bio,scale="percentile_unit",beta=rr["coef"],se=rr["se(coef)"],HR_percentile_unit=rr["exp(coef)"],HR_IQR=exp(rr["coef"]*.5),p_value=rr["Pr(>|z|)"],ph_score_p=cz$table[sc,"p"],ph_global_p=cz$table["GLOBAL","p"],n=nrow(dd),events=sum(dd$event),role="continuity_sensitivity")
    }
  }
}

rmst <- do.call(rbind,rmst_rows); fixed <- do.call(rbind,fixed_rows); tv <- do.call(rbind,tv_rows); cf <- do.call(rbind,coef_rows); qa <- do.call(rbind,qa_rows); cox <- do.call(rbind,cox_rows)
writet(rmst,"04_non_ph_survival/RMST_CONTEXT_MODELS.tsv")
writet(fixed,"04_non_ph_survival/FIXED_TIME_SURVIVAL_MODELS.tsv")
writet(cf,"04_non_ph_survival/FLEXIBLE_PARAMETRIC_MODELS.tsv")
writet(tv,"04_non_ph_survival/TIME_VARYING_INTERACTION.tsv")
writet(qa,"04_non_ph_survival/NON_PH_MODEL_QA.tsv")
writet(cox,"04_non_ph_survival/COX_CONTINUITY_MODELS.tsv")
writet(rmst,"FINAL_RMST_CONTEXT_COMPARISON_V1_7.tsv")
writet(tv,"FINAL_TIME_VARYING_CONTEXT_MODELS_V1_7.tsv")
writeLines(capture.output(sessionInfo()),file.path(run,"session_info/NON_PH_SESSION_INFO.txt"))
cat("completed",nrow(rmst),"RMST rows\n")
}

if(Sys.getenv("NONPH_MODE","MAIN")=="MAIN") run_main()
