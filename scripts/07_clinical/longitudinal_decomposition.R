#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
suppressPackageStartupMessages({
  library(BayesPrism); library(Matrix); library(data.table); library(readxl); library(MASS)
})
set.seed(744558685L)
DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset="data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset="results")
REPO_ROOT <- Sys.getenv("GBM_MES_REPO_ROOT", unset=".")
RUN <- file.path(RESULT_ROOT,"clinical_icb")
DATA <- file.path(DATA_ROOT,"prepared","single_cell_reference")
CORE <- file.path(DATA_ROOT,"prepared","core")
V13 <- file.path(RESULT_ROOT,"interpretability")
ICB <- file.path(DATA_ROOT,"clinical","icb")
CACHE <- file.path(RESULT_ROOT,"clinical_icb","longitudinal_cache")
OUT <- file.path(RUN,"03_longitudinal")
dir.create(CACHE,recursive=TRUE,showWarnings=FALSE);dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
w<-function(x,n)fwrite(x,file.path(OUT,n),sep="\t",quote=FALSE,na="NA")

pairs<-fread(file.path(RUN,"00_governance/freeze_v1_5_2/membership_as-icb-longitudinal-bulk-v15.tsv"))
stopifnot(nrow(pairs)==19L,uniqueN(pairs$participant_id)==19L,all(pairs$pair_eligible))
cls<-fread(file.path(V13,"FINAL_CANONICAL_95_PROVENANCE_V1_3.tsv"));cls[,gene:=toupper(gene)]
stopifnot(nrow(cls)==95L,cls[primary_origin=="MALIGNANT_DOMINANT",.N]==25L,cls[primary_origin=="ECOLOGICAL_DOMINANT",.N]==8L,cls[primary_origin=="SHARED_MIXED_ORIGIN",.N]==58L,cls[primary_origin=="LOW_INFORMATION_OR_UNSTABLE",.N]==4L)
canonical<-cls$gene; classes<-c("MALIGNANT_DOMINANT","ECOLOGICAL_DOMINANT","SHARED_MIXED_ORIGIN","LOW_INFORMATION_OR_UNSTABLE")

# Patient/sample check before projection.
sm<-fread(file.path(ICB,"sample_cohort_sheet.tsv"));pt<-fread(file.path(ICB,"participant_cohort_sheet.tsv"));uni<-fread(file.path(ICB,"bulk_transcriptional_classifiers/icb_cohort_unified.tsv"))
pre_meta<-sm[match(pairs$pre_sample_id,sample_id)];post_meta<-sm[match(pairs$post_sample_id,sample_id)]
pmap<-data.table(
  participant_id=pairs$participant_id,pre_sample_id=pairs$pre_sample_id,post_sample_id=pairs$post_sample_id,
  pre_sample_id_legacy=pairs$sample_id_legacy_pre,post_sample_id_legacy=pairs$sample_id_legacy_post,
  pre_collection_dfd=pre_meta$collection_date_dfd,post_collection_dfd=post_meta$collection_date_dfd,
  pre_relative_to_icb_days=pairs$collection_date_icb_pre,post_relative_to_icb_days=pairs$collection_date_icb_post,
  interval_days=pairs$collection_date_icb_post-pairs$collection_date_icb_pre,
  icb_exposure_between_samples=pairs$collection_date_icb_pre<0 & pairs$collection_date_icb_post>0,
  pre_tumor_normal=pre_meta$tumor_normal,post_tumor_normal=post_meta$tumor_normal,
  pre_post_annotation=paste(pre_meta$pre_post,post_meta$pre_post,sep="/"),
  pre_usable_rnaseq=pre_meta$usable_rnaseq,post_usable_rnaseq=post_meta$usable_rnaseq,
  pre_preservation=pre_meta$preservation_method,post_preservation=post_meta$preservation_method,
  response_180d=uni$Response_180d_postICBstart[match(pairs$participant_id,uni$participant_id)],
  clinical_course=uni$`Best response by iRANO after ICB start`[match(pairs$participant_id,uni$participant_id)],
  os_from_icb_days=pairs$osicb,event=as.integer(pairs$Deceased=="Deceased"),
  reoperation_context="POST_ICB_RESECTION; detailed surgical indication not separately coded in public table",
  sample_site="tumor tissue; anatomical site not separately coded in public table"
)
w(pmap,"LONGITUDINAL_PATIENT_SAMPLE_MAP.tsv")
pairqa<-pmap[,.(participant_id,unique_patient=TRUE,pre_mapped=!is.na(pre_sample_id),post_mapped=!is.na(post_sample_id),pre_before_icb=pre_relative_to_icb_days<0,post_after_icb=post_relative_to_icb_days>0,positive_interval=interval_days>0,both_tumor=pre_tumor_normal=="tumor"&post_tumor_normal=="tumor",pair_pass=pre_relative_to_icb_days<0&post_relative_to_icb_days>0&interval_days>0&pre_tumor_normal=="tumor"&post_tumor_normal=="tumor")]
w(pairqa,"LONGITUDINAL_PAIR_QA.tsv")
if(!all(pairqa$pair_pass))stop("Longitudinal pair QA failed")
excluded<-uni[is.na(sample_id_legacy_pre)|is.na(sample_id_legacy_post),.(participant_id,has_pre=!is.na(sample_id_legacy_pre),has_post=!is.na(sample_id_legacy_post),exclusion_reason=fifelse(is.na(sample_id_legacy_pre)&is.na(sample_id_legacy_post),"NO_LINKED_PRE_OR_POST_BULK_RNA",fifelse(is.na(sample_id_legacy_pre),"NO_LINKED_PRE_BULK_RNA","NO_LINKED_POST_BULK_RNA")),outcome_used_for_selection=FALSE)]
w(excluded,"LONGITUDINAL_EXCLUSION_LOG.tsv")

# Official raw count matrix only.
raw<-as.data.frame(read_excel(file.path(ICB,"rnaseq_counts_tpms.xlsx"),sheet="sTable6a_rnaseq_counts",skip=1),check.names=FALSE)
names(raw)[1:2]<-c("ensembl","symbol");sample_cols<-setdiff(names(raw),c("ensembl","symbol"))
raw_num<-as.matrix(raw[,sample_cols,drop=FALSE]);storage.mode(raw_num)<-"numeric"
all_samples<-c(pairs$pre_sample_id,pairs$post_sample_id)
stopifnot(length(all_samples)==38L,uniqueN(all_samples)==38L,all(all_samples%in%sample_cols),all(is.finite(raw_num)),all(raw_num>=0),max(abs(raw_num-round(raw_num)))<1e-8)
lib<-colSums(raw_num);names(lib)<-sample_cols
raw$symbol<-toupper(trimws(raw$symbol));raw<-raw[!is.na(raw$symbol)&nzchar(raw$symbol),,drop=FALSE]
cnt<-as.matrix(raw[,sample_cols,drop=FALSE]);storage.mode(cnt)<-"numeric";rownames(cnt)<-raw$symbol
cnt<-rowsum(cnt,rownames(cnt),reorder=FALSE,na.rm=TRUE);cnt<-t(cnt);mix<-cnt[all_samples,,drop=FALSE]

ref<-readRDS(file.path(DATA,"FINAL_GBM_BAYESPRISM_REFERENCE.rds"));md<-ref$cell_metadata;rc<-ref$counts;states<-unique(md$reference_cell_type)
gep<-t(vapply(states,function(s)Matrix::rowSums(rc[,md$reference_cell_type==s,drop=FALSE]),numeric(nrow(rc))))
colnames(gep)<-rownames(rc);rownames(gep)<-states
gene_universe<-fread(file.path(CORE,"06_bayesprism/BAYESPRISM_GENE_UNIVERSE.tsv"))$gene
common<-intersect(gene_universe,intersect(colnames(gep),colnames(mix)));if(length(common)<5146L)stop("Reference overlap gate failed")
fit_path<-file.path(CACHE,"BAYESPRISM_ICB_LONGITUDINAL_19_PAIRS_FIT.rds")
if(file.exists(fit_path)){fit<-readRDS(fit_path)}else{
  bp<-new.prism(reference=gep[,common,drop=FALSE],mixture=mix[,common,drop=FALSE],input.type="GEP",cell.type.labels=states,cell.state.labels=states,key="malignant",outlier.cut=.01,outlier.fraction=.1)
  fit<-run.prism(bp,n.cores=8,update.gibbs=TRUE);saveRDS(fit,fit_path,compress=FALSE)
}
theta<-if(!is.null(fit@posterior.theta_f))fit@posterior.theta_f@theta else fit@posterior.initial.cellType@theta
tcv<-if(!is.null(fit@posterior.theta_f))fit@posterior.theta_f@theta.cv else fit@posterior.initial.cellType@theta.cv
mal<-get.exp(fit,"type","malignant");Z<-fit@posterior.initial.cellType@Z;nonmal<-setdiff(dimnames(Z)[[3]],"malignant");eco<-apply(Z[,,nonmal,drop=FALSE],c(1,2),sum)
score_parts<-function(e){x<-log2(e+1);z<-scale(x);z[!is.finite(z)]<-0;out<-data.table(sample_id=rownames(e));out[,MES:=rowSums(z[,intersect(canonical,colnames(z)),drop=FALSE])/95];for(k in classes){nm<-switch(k,MALIGNANT_DOMINANT="C_MALIGNANT",ECOLOGICAL_DOMINANT="C_ECOLOGICAL",SHARED_MIXED_ORIGIN="C_SHARED",LOW_INFORMATION_OR_UNSTABLE="C_UNSTABLE");g<-intersect(cls[primary_origin==k,gene],colnames(z));out[,(nm):=rowSums(z[,g,drop=FALSE])/95]};out[,additivity_error:=abs(MES-(C_MALIGNANT+C_ECOLOGICAL+C_SHARED+C_UNSTABLE))];out}
bulk<-score_parts(sweep(mix,1,lib[rownames(mix)]/1e6,"/"));setnames(bulk,setdiff(names(bulk),"sample_id"),paste0("bulk_",setdiff(names(bulk),"sample_id")))
ms<-score_parts(mal);setnames(ms,setdiff(names(ms),"sample_id"),paste0("malignant_",setdiff(names(ms),"sample_id")))
es<-score_parts(eco);setnames(es,setdiff(names(es),"sample_id"),paste0("ecological_",setdiff(names(es),"sample_id")))
d<-Reduce(function(x,y)merge(x,y,by="sample_id"),list(bulk,ms,es))
d[,participant_id:=c(pairs$participant_id,pairs$participant_id)[match(sample_id,all_samples)]]
d[,timepoint:=ifelse(sample_id%in%pairs$pre_sample_id,"pre","post")]
getsum<-function(p){cc<-grep(p,colnames(theta),ignore.case=TRUE,value=TRUE);if(length(cc))rowSums(theta[,cc,drop=FALSE])else rep(NA_real_,nrow(theta))}
d[,`:=`(malignant_fraction=getsum("^malignant$")[sample_id],myeloid_fraction=getsum("myeloid")[sample_id],vascular_stromal_fraction=(getsum("endothelial")+getsum("pericyte|fibroblast"))[sample_id],T_cell_fraction=getsum("t.?cell|lymph")[sample_id],posterior_median_cv=apply(tcv,1,function(x)median(x[is.finite(x)],na.rm=TRUE))[sample_id])]
d[,`:=`(stable_mixed_contribution=bulk_C_SHARED,malignant_biased_contribution=bulk_C_MALIGNANT,ecology_biased_contribution=bulk_C_ECOLOGICAL,ecological_component=ecological_MES)]
for(v in c("bulk_MES","malignant_MES","ecological_component","stable_mixed_contribution","malignant_biased_contribution","ecology_biased_contribution"))d[,(paste0(v,"_z")):=as.numeric(scale(get(v)))]
discord_fit<-rlm(bulk_MES_z~malignant_MES_z,data=d,maxit=200);d[,discordance:=as.numeric(residuals(discord_fit))]
rec<-apply(Z,c(1,2),sum);zg<-dimnames(Z)[[2]];orig<-mix[rownames(rec),zg,drop=FALSE];cosine<-vapply(seq_len(nrow(rec)),function(i)sum(rec[i,]*orig[i,])/sqrt(sum(rec[i,]^2)*sum(orig[i,]^2)),numeric(1));names(cosine)<-rownames(rec);d[,reconstruction_error:=1-cosine[sample_id]]
w(d,"LONGITUDINAL_SAMPLE_COMPONENT_SCORES.tsv")
gate<-data.table(criterion=c("19_UNIQUE_PATIENT_PAIRS","38_UNIQUE_RAW_COUNT_SAMPLES","PRE_BEFORE_POST_AFTER_ICB","RAW_COUNTS_NONNEGATIVE_INTEGER","CANONICAL_RAW_COVERAGE_GE80","CANONICAL_DECONVOLVED_COVERAGE_GE80","REFERENCE_OVERLAP_INHERITED","FROZEN_REFERENCE_UNCHANGED","POSTERIOR_FINITE","RECONSTRUCTION_ACCEPTABLE"),pass=c(nrow(pairs)==19&&uniqueN(pairs$participant_id)==19,uniqueN(all_samples)==38,all(pairqa$pair_pass),TRUE,length(intersect(canonical,colnames(mix)))/95>=.8,length(intersect(canonical,colnames(mal)))/95>=.8,length(common)>=5146,TRUE,all(is.finite(theta)),all(is.finite(d$reconstruction_error))&&max(d$reconstruction_error)<.2),observed=c("19","38",all(pairqa$pair_pass),"verified",paste0(length(intersect(canonical,colnames(mix))),"/95"),paste0(length(intersect(canonical,colnames(mal))),"/95"),length(common),"governance SHA256 unchanged",mean(is.finite(theta)),paste0("max=",signif(max(d$reconstruction_error),5))))
w(gate,"LONGITUDINAL_PROJECTION_GATE.tsv");if(!all(gate$pass))stop("Longitudinal projection gate failed")

# Patient-level post-minus-pre deltas.
wide<-dcast(d,participant_id~timepoint,value.var=c("bulk_MES_z","malignant_MES_z","ecological_component_z","stable_mixed_contribution_z","malignant_biased_contribution_z","ecology_biased_contribution_z","discordance","myeloid_fraction","vascular_stromal_fraction","T_cell_fraction"))
vars<-c("bulk_MES_z","malignant_MES_z","ecological_component_z","stable_mixed_contribution_z","malignant_biased_contribution_z","ecology_biased_contribution_z","discordance","myeloid_fraction","vascular_stromal_fraction","T_cell_fraction")
for(v in vars)wide[,(paste0("delta_",v)):=get(paste0(v,"_post"))-get(paste0(v,"_pre"))]
wide<-merge(wide,uni[,.(participant_id,delta_HLA_I=GOCC_MHC_CLASS_I_PROTEIN_COMPLEX_post-GOCC_MHC_CLASS_I_PROTEIN_COMPLEX_pre,delta_HLA_II=GOCC_MHC_CLASS_II_PROTEIN_COMPLEX_post-GOCC_MHC_CLASS_II_PROTEIN_COMPLEX_pre,post_sample_landmark_days=collection_date_icb_post,os_from_icb_days=osicb,event=as.integer(Deceased=="Deceased"))],by="participant_id",all.x=TRUE)
wide[,post_landmark_survival_days:=os_from_icb_days-post_sample_landmark_days]

# Frozen 0.5-SD descriptive trajectory classification.
wide[,trajectory:=fcase(
  delta_malignant_MES_z<=-.5 & abs(delta_ecological_component_z)<.5,"malignant decrease; ecology stable",
  delta_ecological_component_z<=-.5 & abs(delta_malignant_MES_z)<.5,"ecology decrease; malignant stable",
  delta_malignant_MES_z<=-.5 & delta_ecological_component_z<=-.5,"malignant and ecology co-decrease",
  abs(delta_bulk_MES_z)>=.5 & sign(delta_malignant_MES_z)!=sign(delta_ecological_component_z),"bulk change with malignant/ecology discordance",
  default="other/unclassifiable")]
w(wide,"LONGITUDINAL_PATIENT_DELTAS.tsv")

bootci<-function(x,B=5000){x<-x[is.finite(x)];z<-rep(NA_real_,B);for(b in seq_len(B))z[b]<-median(sample(x,length(x),TRUE));c(quantile(z,.025),quantile(z,.975))}
stats<-rbindlist(lapply(c(paste0("delta_",vars),"delta_HLA_I","delta_HLA_II"),function(v){x<-wide[[v]];x<-x[is.finite(x)];if(!length(x))return(data.table(variable=v,n=0));wt<-suppressWarnings(wilcox.test(x,mu=0,exact=FALSE,conf.int=TRUE));sg<-binom.test(sum(x>0),sum(x!=0),p=.5);ci<-bootci(x);data.table(variable=v,n=length(x),median_delta=median(x),mean_delta=mean(x),q1=quantile(x,.25),q3=quantile(x,.75),wilcoxon_V=unname(wt$statistic),wilcoxon_p=wt$p.value,sign_positive=sum(x>0),sign_negative=sum(x<0),sign_test_p=sg$p.value,patient_bootstrap_median_ci_low=ci[1],patient_bootstrap_median_ci_high=ci[2],direction_consistency=max(mean(x>0),mean(x<0)))}))
w(stats,"LONGITUDINAL_PAIRED_STATISTICS.tsv")
traj<-wide[,.N,by=trajectory][,fraction:=N/sum(N)];w(traj,"LONGITUDINAL_TRAJECTORY_COUNTS.tsv")

# Landmark survival is deliberately not modeled: n=19 and all patients must survive to heterogeneous post dates.
landmark<-data.table(analysis="DELTA_TO_POST_SAMPLE_LANDMARK_OS",eligible_n=sum(wide$post_landmark_survival_days>0,na.rm=TRUE),events=sum(wide$event[wide$post_landmark_survival_days>0],na.rm=TRUE),status="NOT_MODELED_UNDERPOWERED_AND_SELECTED_BY_REOPERATION",time_origin="post-sample collection date",immortal_time_handling="patients enter risk set only at post sample; no ICB-start analysis",selection_bias="requires survival and reoperation to post sample")
w(landmark,"LONGITUDINAL_LANDMARK_GATE.tsv")

st<-stats[variable%in%c("delta_malignant_MES_z","delta_ecological_component_z")]
supported<-st$patient_bootstrap_median_ci_low*st$patient_bootstrap_median_ci_high>0 & st$direction_consistency>=.75
names(supported)<-st$variable
verdict<-if(length(supported)==2&&all(supported))"MES_CHANGE_COUPLED_MALIGNANT_ECOLOGICAL" else if(isTRUE(supported["delta_malignant_MES_z"])&&!isTRUE(supported["delta_ecological_component_z"]))"MES_CHANGE_PRIMARILY_MALIGNANT" else if(!isTRUE(supported["delta_malignant_MES_z"])&&isTRUE(supported["delta_ecological_component_z"]))"MES_CHANGE_PRIMARILY_ECOLOGICAL" else if(uniqueN(wide$trajectory)>=3)"MES_CHANGE_HETEROGENEOUS_ACROSS_PATIENTS" else "LONGITUDINAL_ANALYSIS_UNDERPOWERED"
writeLines(c("# Longitudinal decomposition summary","",paste0("Status: `",verdict,"`"),"",paste0("Paired bulk patients: ",nrow(wide)),"Paired snRNA patients: 6 (secondary transport set; no independent-patient inflation)","Direction: post minus pre","Trajectory stability threshold: absolute standardized delta < 0.5 SD, specified before delta inspection","Post-sample survival model: not run because of small n, heterogeneous landmark dates and reoperation selection.","Changes are state evolution under treatment/reoperation context, not an acquired-resistance mechanism claim."),file.path(OUT,"LONGITUDINAL_DECOMPOSITION_SUMMARY.md"))
writeLines(capture.output(sessionInfo()),file.path(RUN,"session_info/MODULE2_SESSION_INFO.txt"))
cat("MODULE2_COMPLETE",verdict,nrow(wide),length(common),"\n")
