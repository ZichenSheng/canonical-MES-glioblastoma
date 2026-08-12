#!/usr/bin/env Rscript
Sys.setenv(NONPH_MODE="LIB")
source(file.path(Sys.getenv("GBM_MES_REPO_ROOT", unset="."),"scripts","07_clinical","nonproportional_hazards.R"))
set.seed(17032026L)
times_grid <- c(1,3,6,9,12,18,24)

sens <- list(); concord <- list(); anti_rmst <- list(); anti_fixed <- list(); anti_tv <- list()
for(stage0 in c("newly_diagnosed","recurrent")) {
  tau <- if(stage0=="newly_diagnosed") 18 else 12
  di0 <- dall[dall$stage==stage0 & dall$context=="ICB",,drop=FALSE]
  dg <- dall[dall$stage==stage0 & dall$context=="GLASS",,drop=FALSE]
  di <- di0[grepl("nivolumab|pembrolizumab",tolower(di0$agent)),,drop=FALSE]
  for(bio in c("stable_mixed","bulk_MES")) {
    sc <- paste0(bio,"_percentile")
    fi <- fit_one(di,sc); fg <- fit_one(dg,sc)
    if(!is.null(fi) && !is.null(fg)) {
      ei <- effects(fi,tau); eg <- effects(fg,tau); boot <- boot_pair(di,dg,sc,tau,times_grid)
      ci <- qci(boot[,"rmst_delta"])
      anti_rmst[[length(anti_rmst)+1]] <- data.frame(sensitivity="anti_PD1_only",stage=stage0,biomarker=bio,tau_months=tau,icb_n=nrow(di),icb_events=sum(di$event),glass_n=nrow(dg),glass_events=sum(dg$event),icb_effect=ei$rmst,glass_effect=eg$rmst,context_difference=ei$rmst-eg$rmst,ci_low=ci[1],ci_high=ci[2],p_value=ci[3],bootstrap_success=ci[4])
      for(j in 1:2) {
        tt<-c(12,18)[j]; cj<-qci(boot[,c("fixed12_delta","fixed18_delta")[j]])
        anti_fixed[[length(anti_fixed)+1]] <- data.frame(sensitivity="anti_PD1_only",stage=stage0,biomarker=bio,time_months=tt,icb_effect=ei$fixed[j],glass_effect=eg$fixed[j],context_difference=ei$fixed[j]-eg$fixed[j],ci_low=cj[1],ci_high=cj[2],p_value=cj[3],icb_at_risk=sum(di$time_months>=tt),glass_at_risk=sum(dg$time_months>=tt),support=ifelse(sum(di$time_months>=tt)>=5 & sum(dg$time_months>=tt)>=5,"ADEQUATE","LIMITED"))
      }
      li<-loghr_time(fi,times_grid); lg<-loghr_time(fg,times_grid)
      for(j in seq_along(times_grid)) {cj<-qci(boot[,paste0("tv",times_grid[j])]); anti_tv[[length(anti_tv)+1]]<-data.frame(sensitivity="anti_PD1_only",stage=stage0,biomarker=bio,time_months=times_grid[j],delta_logHR=li[j]-lg[j],ci_low=cj[1],ci_high=cj[2],p_value=cj[3],bootstrap_success=cj[4])}
    }
    for(ctx in c("ICB","GLASS")) {
      dd <- if(ctx=="ICB") di0 else dg
      for(scale0 in c("percentile","rank_normal","sd")) {
        scol <- paste0(bio,"_",scale0)
        dd$age_z <- if(sd(dd$age,na.rm=TRUE)>0) (dd$age-mean(dd$age,na.rm=TRUE))/sd(dd$age,na.rm=TRUE) else 0
        fit <- coxph(as.formula(paste0("Surv(time_months,event)~",scol,"+age_z")),data=dd,x=TRUE)
        z<-cox.zph(fit); rr<-summary(fit)$coef[scol,]
        contrast <- if(scale0=="percentile") .5 else if(scale0=="rank_normal") qnorm(.75)-qnorm(.25) else 1
        sens[[length(sens)+1]] <- data.frame(module="score_scale",stage=stage0,context=ctx,biomarker=bio,sensitivity=scale0,estimate=rr["coef"]*contrast,effect_measure="logHR_contrast",HR=exp(rr["coef"]*contrast),ci_low=exp((rr["coef"]-1.96*rr["se(coef)"])*contrast),ci_high=exp((rr["coef"]+1.96*rr["se(coef)"])*contrast),p_value=rr["Pr(>|z|)"],ph_p=z$table[scol,"p"],n=nrow(dd),events=sum(dd$event),status="ESTIMATED")
      }
    }
  }
}

ar<-do.call(rbind,anti_rmst); af<-do.call(rbind,anti_fixed); atv<-do.call(rbind,anti_tv); ss<-do.call(rbind,sens)
writet(ar,"06_sensitivity/ANTI_PD1_RMST.tsv"); writet(af,"06_sensitivity/ANTI_PD1_FIXED_TIME.tsv"); writet(atv,"06_sensitivity/ANTI_PD1_FLEXIBLE_TIME_VARYING.tsv")

tz <- read.delim(file.path(run,"02_time_zero_alignment/SAMPLE_TO_TREATMENT_INTERVAL.tsv"),check.names=FALSE)
extra <- list()
for(stage0 in c("newly_diagnosed","recurrent")) for(w in c(28,42,60)) {
  ii<-tz[tz$dataset=="ICB" & tz$stage==stage0,]; gg<-tz[tz$dataset=="GLASS" & tz$stage==stage0,]
  extra[[length(extra)+1]]<-data.frame(module="sample_window",stage=stage0,context="ICB_vs_GLASS",biomarker="stable_mixed",sensitivity=paste0("0_",w,"_days"),estimate=NA,effect_measure="context_difference",HR=NA,ci_low=NA,ci_high=NA,p_value=NA,ph_p=NA,n=sum(ii[[paste0("window_0_",w)]]=="YES"),events=NA,status="NOT_EVALUABLE_GLASS_EXACT_TREATMENT_START_MISSING")
}
extra[[length(extra)+1]]<-data.frame(module="weighting",stage="newly_diagnosed",context="ICB_vs_GLASS",biomarker="stable_mixed",sensitivity="overlap_vs_IPTW_vs_unweighted",estimate=NA,effect_measure="interaction",HR=NA,ci_low=NA,ci_high=NA,p_value=NA,ph_p=NA,n=NA,events=NA,status="NOT_APPLICABLE_TIME_ZERO_GATE")
extra[[length(extra)+1]]<-data.frame(module="weighting",stage="recurrent",context="ICB_vs_GLASS",biomarker="stable_mixed",sensitivity="overlap_vs_IPTW_vs_unweighted",estimate=NA,effect_measure="interaction",HR=NA,ci_low=NA,ci_high=NA,p_value=NA,ph_p=NA,n=NA,events=NA,status="NOT_APPLICABLE_TIME_ZERO_GATE")
extra[[length(extra)+1]]<-data.frame(module="leave_one_source_out",stage="all",context="ICB",biomarker="stable_mixed",sensitivity="source",estimate=NA,effect_measure="interaction",HR=NA,ci_low=NA,ci_high=NA,p_value=NA,ph_p=NA,n=40,events=38,status="NOT_EVALUABLE_SINGLE_ICB_SOURCE_LABEL")
extra[[length(extra)+1]]<-data.frame(module="calendar_era",stage="all",context="GLASS",biomarker="stable_mixed",sensitivity="prespecified_era",estimate=NA,effect_measure="interaction",HR=NA,ci_low=NA,ci_high=NA,p_value=NA,ph_p=NA,n=237,events=215,status="NOT_EVALUABLE_ABSOLUTE_DATES_MISSING")
all <- rbind(ss,do.call(rbind,extra))
writet(all,"06_sensitivity/PRESPECIFIED_SENSITIVITY_RESULTS.tsv")

cm <- data.frame(check=c("full_ICB_newly_RMST","anti_PD1_newly_RMST","full_ICB_recurrent_RMST","anti_PD1_recurrent_RMST","score_scale_direction","weighting_concordance","time_zero","overall"),result=c("directionally_positive_CI_crosses_zero","directionally_positive_CI_crosses_zero","imprecise","imprecise","generally_directionally_consistent_but_PH_sensitive","not_evaluable","not_alignable","NON_PH_CONTEXT_SIGNAL_IMPRECISE"),claim_impact=c("no upgrade","no upgrade","no upgrade","no upgrade","supports continuity only","target trial stopped","causal claim prohibited","retain observational boundary"),stringsAsFactors=FALSE)
writet(cm,"06_sensitivity/SENSITIVITY_CONCORDANCE_MATRIX.tsv"); writet(cm,"FINAL_SENSITIVITY_MATRIX_V1_7.tsv")
writeLines(c("# Sensitivity gate","","`NON_PH_CONTEXT_SIGNAL_IMPRECISE`","","The newly diagnosed direction is similar in the full ICB and anti-PD-1-only analyses, but interval estimates and time-varying coefficients are unstable. Recurrent analyses are imprecise. Window and weighting sensitivities are not evaluable because GLASS lacks authenticated treatment start. No additional subgroup, date window, or propensity algorithm was introduced."),file.path(run,"06_sensitivity/SENSITIVITY_GATE.md"))
writeLines(capture.output(sessionInfo()),file.path(run,"session_info/SENSITIVITY_SESSION_INFO.txt"))
cat("sensitivity complete\n")
