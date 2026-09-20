#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(MASS)})
set.seed(744558686L)
DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset="data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset="results")
REPO_ROOT <- Sys.getenv("GBM_MES_REPO_ROOT", unset=".")
OUT<-file.path(RESULT_ROOT,"clinical_icb","orthogonal");dir.create(OUT,recursive=TRUE,showWarnings=FALSE);w<-function(x,n)fwrite(x,file.path(OUT,n),sep="\t",quote=FALSE,na="NA")
d<-fread(file.path(OUT,"HLA_CD3_PATIENT_LEVEL_TABLE.tsv"))

spec<-data.table(
  endpoint_id=c("O1","O2","O3","O4","O5","S1","S2","S3"),
  tier=c(rep("PRIMARY",5),rep("SECONDARY",3)),
  x=c("malignant_MES","ecological_component","stable_mixed_contribution","stable_mixed_contribution","discordance","malignant_MES","myeloid_fraction","tumor_cell_HLA_I"),
  y=c("tumor_cell_HLA_I","cd3_density","tumor_cell_HLA_I","cd3_density","cd3_density","tumor_cell_HLA_II","cd3_density","tumor_cell_HLA_II"),
  contrast=c("malignant MES vs tumour-cell HLA-I","ecological component vs CD3 density","stable mixed vs tumour-cell HLA-I","stable mixed vs CD3 density","discordance vs CD3 density","malignant MES vs tumour-cell HLA-II","myeloid fraction vs CD3 density","HLA-I vs HLA-II tumour phenotype")
)
bootrho<-function(x,y,B=5000){ok<-is.finite(x)&is.finite(y);x<-x[ok];y<-y[ok];n<-length(x);if(n<5)return(c(NA,NA,NA));rho<-cor(x,y,method="spearman");z<-rep(NA_real_,B);for(b in seq_len(B)){ii<-sample.int(n,n,TRUE);z[b]<-suppressWarnings(cor(x[ii],y[ii],method="spearman"))};c(rho,quantile(z,.025,na.rm=TRUE),quantile(z,.975,na.rm=TRUE))}
rows<-list()
for(i in seq_len(nrow(spec))){s<-spec[i];x<-d[[s$x]];y<-d[[s$y]];ok<-is.finite(x)&is.finite(y);n<-sum(ok);r<-bootrho(x,y);ct<-if(n>=3)cor.test(x[ok],y[ok],method="spearman",exact=FALSE)else NULL;rr<-if(n>=8)tryCatch(rlm(scale(y[ok])~scale(x[ok]),maxit=200),error=function(e)NULL)else NULL;rows[[i]]<-data.table(endpoint_id=s$endpoint_id,tier=s$tier,contrast=s$contrast,x=s$x,y=s$y,n=n,spearman_rho=r[1],spearman_p=if(is.null(ct))NA_real_ else ct$p.value,bootstrap_ci_low=r[2],bootstrap_ci_high=r[3],robust_beta=if(is.null(rr))NA_real_ else coef(rr)[2],robust_se=if(is.null(rr))NA_real_ else summary(rr)$coefficients[2,2],robust_p=if(is.null(rr))NA_real_ else 2*pnorm(-abs(coef(rr)[2]/summary(rr)$coefficients[2,2])))}
res<-rbindlist(rows,fill=TRUE);res[,BH_FDR:=p.adjust(spearman_p,method="BH")];w(res,"HLA_CD3_ORTHOGONAL_ASSOCIATIONS.tsv")

# Formally compare dependent correlations using patient bootstrap; no significance-vs-nonsignificance shortcut.
depboot<-function(a,b,y,B=5000){ok<-is.finite(a)&is.finite(b)&is.finite(y);a<-a[ok];b<-b[ok];y<-y[ok];n<-length(y);obs<-cor(a,y,method="spearman")-cor(b,y,method="spearman");z<-rep(NA_real_,B);for(k in seq_len(B)){ii<-sample.int(n,n,TRUE);z[k]<-suppressWarnings(cor(a[ii],y[ii],method="spearman")-cor(b[ii],y[ii],method="spearman"))};data.table(n=n,rho_difference=obs,ci_low=quantile(z,.025,na.rm=TRUE),ci_high=quantile(z,.975,na.rm=TRUE),bootstrap_two_sided_p=2*min(mean(z<=0,na.rm=TRUE),mean(z>=0,na.rm=TRUE)))}
cmp1<-depboot(d$stable_mixed_contribution,d$malignant_MES,d$cd3_density);cmp1[,contrast:="rho(shared, CD3) - rho(malignant, CD3)"]
cmp2<-depboot(d$stable_mixed_contribution,d$malignant_MES,d$tumor_cell_HLA_I);cmp2[,contrast:="rho(shared, tumour HLA-I) - rho(malignant, tumour HLA-I)"]
w(rbindlist(list(cmp1,cmp2),fill=TRUE),"HLA_CD3_DEPENDENT_CORRELATION_BOOTSTRAP.tsv")

primary<-res[tier=="PRIMARY"]
supported<-primary[is.finite(bootstrap_ci_low)&bootstrap_ci_low*bootstrap_ci_high>0&BH_FDR<.1]
mal_anchor<-"O1"%in%supported$endpoint_id;eco_anchor<-"O2"%in%supported$endpoint_id;shared_anchor<-any(c("O3","O4")%in%supported$endpoint_id)
verdict<-if(mal_anchor)"MALIGNANT_MES_ORTHOGONALLY_ANCHORED_TO_HLAI" else if(eco_anchor)"ECOLOGICAL_MES_ORTHOGONALLY_ANCHORED_TO_CD3" else if(shared_anchor)"SHARED_MES_ORTHOGONALLY_ANCHORED" else if(any(is.finite(primary$spearman_rho)&abs(primary$spearman_rho)>=.3))"PARTIAL_ORTHOGONAL_ANCHORING" else if(all(primary$n<5))"NOT_EVALUABLE" else "NO_CLEAR_ORTHOGONAL_ANCHOR"
writeLines(c("# HLA/CD3 orthogonal validation summary","",paste0("Status: `",verdict,"`"),"",paste0("Patient-level mIF/IHC rows: ",nrow(d)),paste0("Patients with exact RNA-component tissue match: ",sum(d$exact_rna_component_match,na.rm=TRUE)),"Cell and ROI rows were aggregated before inference; all primary tests use one row per patient.","BH-FDR is applied within this module. Correlation differences use dependent-correlation patient bootstrap.","Spatial nearest-neighbor/contact analysis is not evaluable because supplied raw tables contain ROI labels but no cell x/y coordinates."),file.path(OUT,"HLA_CD3_VALIDATION_SUMMARY.md"))
writeLines(capture.output(sessionInfo()),file.path(RUN,"session_info/MODULE3_SESSION_INFO.txt"))
cat("MODULE3_COMPLETE",verdict,nrow(d),"\n")
