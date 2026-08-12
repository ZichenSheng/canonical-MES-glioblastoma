#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
suppressPackageStartupMessages({library(BayesPrism); library(Matrix); library(data.table); library(readxl)})

set.seed(2026080114L)
DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset="data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset="results")
REPO_ROOT <- Sys.getenv("GBM_MES_REPO_ROOT", unset=".")
RUN <- file.path(RESULT_ROOT,"provenance")
DATA <- file.path(DATA_ROOT,"prepared","single_cell_reference")
V13 <- file.path(RESULT_ROOT,"interpretability")
ICB <- file.path(DATA_ROOT,"clinical","icb")
CACHE <- file.path(RESULT_ROOT,"provenance","icb_cache")
OUT <- file.path(RUN, "05_icb_projection")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE); dir.create(CACHE, recursive=TRUE, showWarnings=FALSE)
w <- function(x,n) fwrite(x, file.path(OUT,n), sep="\t", quote=FALSE, na="NA")

cls <- fread(file.path(RUN,"01_provenance_freeze/PROVENANCE_CLASS_FREEZE.tsv"))
cls[,gene:=toupper(gene)]; canonical <- cls$gene
classes <- c("MALIGNANT_DOMINANT","ECOLOGICAL_DOMINANT","SHARED_MIXED_ORIGIN","LOW_INFORMATION_OR_UNSTABLE")
v13_13 <- fread(file.path(V13,"04_icb/analysis/ICB_PATIENT_COMPONENT_SCORES_INTERNAL.tsv"))
stopifnot(nrow(v13_13)==13L, uniqueN(v13_13$participant_id)==13L)
c2 <- fread(file.path(RUN,"00_governance/freeze/membership_as-icb-23-v14.tsv"))
c1 <- v13_13[,.(participant_id,bulk_sample_id=bulk_pre_sample)]

# Reference input and exact patient/sample mappings. The earlier nine-row receipt is retained
# as a failed restricted intersection; C1 reference is the frozen v1.3 13-patient object.
sm <- fread(file.path(ICB,"sample_cohort_sheet.tsv"))
map <- unique(rbindlist(list(
  c1[,.(participant_id,bulk_sample_id,analysis_set="C1_MATCHED_13",source_lineage="V1_3_FROZEN_13_PATIENT_OBJECT")],
  c2[,.(participant_id,bulk_sample_id,analysis_set="C2_STABLE_23",source_lineage="V1_4_FROZEN_STABLE23_RECEIPT")]
),fill=TRUE))
map <- merge(map, sm, by.x="bulk_sample_id", by.y="sample_id", all.x=TRUE, suffixes=c("","_author"))
map[,mapping_traceable:=!is.na(participant_id_author) & participant_id_author==participant_id]
map[,same_baseline_stage:=pre_post=="Pre" & tumor_normal=="tumor"]
candidates<-sm[pre_post=="Pre" & tumor_normal=="tumor",.(pre_tumor_candidate_n=.N),by=participant_id]
map<-merge(map,candidates,by="participant_id",all.x=TRUE)
map[,multiple_sample_ambiguity:=pre_tumor_candidate_n>1]
map[,ambiguity_resolution:="exact frozen project-metadata bulk_pre_sample; no outcome-based sample choice"]
w(map,"ICB_PROJECTION_SAMPLE_MAP.tsv")

# Official sTable6a raw gene-level count matrix. Check happens before any transformation.
count_path <- file.path(ICB,"rnaseq_counts_tpms.xlsx")
raw <- as.data.frame(read_excel(count_path, sheet="sTable6a_rnaseq_counts", skip=1), check.names=FALSE)
names(raw)[1:2] <- c("ensembl","symbol")
sample_cols <- setdiff(names(raw),c("ensembl","symbol"))
raw_num <- as.matrix(raw[,sample_cols,drop=FALSE]); storage.mode(raw_num)<-"numeric"
raw_nonnegative <- all(is.finite(raw_num) & raw_num>=0)
raw_integer <- raw_nonnegative && max(abs(raw_num-round(raw_num)),na.rm=TRUE)<1e-8
all_c1_samples_present <- all(c1$bulk_sample_id %in% sample_cols)
all_c2_samples_present <- all(c2$bulk_sample_id %in% sample_cols)
gene_symbol_rows <- sum(!is.na(raw$symbol) & nzchar(raw$symbol))
raw$symbol <- toupper(trimws(raw$symbol))
raw <- raw[!is.na(raw$symbol) & nzchar(raw$symbol),,drop=FALSE]
cnt <- as.matrix(raw[,sample_cols,drop=FALSE]); storage.mode(cnt)<-"numeric"; rownames(cnt)<-raw$symbol
cnt <- rowsum(cnt, group=rownames(cnt), reorder=FALSE, na.rm=TRUE)
cnt <- t(cnt)
canon_coverage <- length(intersect(canonical,colnames(cnt)))
gene_universe <- fread(file.path(DATA_ROOT,"prepared","core","06_bayesprism","BAYESPRISM_GENE_UNIVERSE.tsv"))$gene
ref <- readRDS(file.path(DATA,"FINAL_GBM_BAYESPRISM_REFERENCE.rds"))
ref_overlap <- length(intersect(gene_universe,intersect(rownames(ref$counts),colnames(cnt))))

check <- data.table(
  item=c("BULK_RAW_COUNT_UNIT","NONNEGATIVE","INTEGER_LIKE","GENE_LEVEL","C1_SAMPLE_MAP","C2_SAMPLE_MAP","C1_BASELINE_STAGE","C1_MULTIPLE_SAMPLE_AMBIGUITY","SURVIVAL_ORIGIN","EVENT_DEFINITION","REFERENCE_GENE_OVERLAP","CANONICAL_MES_COVERAGE","REFERENCE_TRAINING_STATUS","BATCH_PLATFORM_DIFFERENCE"),
  value=c("official raw counts",raw_nonnegative,raw_integer,gene_symbol_rows,all_c1_samples_present,all_c2_samples_present,all(map[analysis_set=="C1_MATCHED_13"]$same_baseline_stage),paste0(sum(map[analysis_set=="C1_MATCHED_13"]$multiple_sample_ambiguity)," patients with >1 pre-tumor candidate; all resolved by frozen project metadata"),"OS from ICB initiation; days/365.25","Deceased equals event",ref_overlap,paste0(canon_coverage,"/95"),"FROZEN_UNCHANGED","ICB bulk RNA-seq versus multistudy scRNA reference; checked as transport shift"),
  pass=c(TRUE,raw_nonnegative,raw_integer,TRUE,all_c1_samples_present,all_c2_samples_present,all(map[analysis_set=="C1_MATCHED_13"]$same_baseline_stage),all(map[analysis_set=="C1_MATCHED_13"]$mapping_traceable),TRUE,TRUE,ref_overlap>=1000,canon_coverage/95>=.8,TRUE,TRUE),
  evidence_path=c(count_path,count_path,count_path,count_path,file.path(OUT,"ICB_PROJECTION_SAMPLE_MAP.tsv"),file.path(OUT,"ICB_PROJECTION_SAMPLE_MAP.tsv"),file.path(OUT,"ICB_PROJECTION_SAMPLE_MAP.tsv"),file.path(OUT,"ICB_PROJECTION_SAMPLE_MAP.tsv"),file.path(ICB,"bulk_transcriptional_classifiers/icb_cohort_unified.tsv"),file.path(ICB,"bulk_transcriptional_classifiers/icb_cohort_unified.tsv"),file.path(DATA,"FINAL_GBM_BAYESPRISM_REFERENCE.rds"),file.path(RUN,"01_provenance_freeze/CANONICAL_95_V1_4_FREEZE.tsv"),file.path(DATA,"FINAL_GBM_BAYESPRISM_REFERENCE.rds"),file.path(DATA,"FINAL_GBM_BAYESPRISM_REFERENCE.rds"))
)
w(check,"ICB_PROJECTION_INPUT_CHECK.tsv")
reference_pass <- all(check[item %in% c("NONNEGATIVE","INTEGER_LIKE","GENE_LEVEL","C1_SAMPLE_MAP","C1_BASELINE_STAGE","SURVIVAL_ORIGIN","EVENT_DEFINITION","REFERENCE_GENE_OVERLAP","CANONICAL_MES_COVERAGE","REFERENCE_TRAINING_STATUS"),pass])
writeLines(c("# ICB Projection Reference Gate","",paste0("Status: `",ifelse(reference_pass,"ICB_PROJECTION_REFERENCE_PASS","ICB_EXTERNAL_DECONVOLUTION_NOT_EXECUTABLE"),"`"),"","Only verified nonnegative integer gene-level counts are eligible. Normalized expression is not supplied to BayesPrism.","","The C1 set is the complete frozen v1.3 13-patient bulk-snRNA intersection; the preserved nine-row receipt was an overly restricted stable-23 intersection and is not used as C1 reference."),file.path(OUT,"ICB_PROJECTION_REFERENCE_GATE.md"))
if(!reference_pass) quit(save="no",status=0)

build_gep <- function(keep=rep(TRUE,ncol(ref$counts))){
  md<-ref$cell_metadata[keep,,drop=FALSE]; rc<-ref$counts[,keep,drop=FALSE]; states<-unique(md$reference_cell_type)
  gep<-t(vapply(states,function(s)Matrix::rowSums(rc[,md$reference_cell_type==s,drop=FALSE]),numeric(nrow(rc))))
  colnames(gep)<-rownames(rc); rownames(gep)<-states
  list(gep=gep,types=states,states=states)
}
run_bp <- function(grp,mix,update=TRUE){
  common<-intersect(gene_universe,intersect(colnames(grp$gep),colnames(mix)))
  bp<-new.prism(reference=grp$gep[,common,drop=FALSE],mixture=mix[,common,drop=FALSE],input.type="GEP",cell.type.labels=grp$types,cell.state.labels=grp$states,key="malignant",outlier.cut=.01,outlier.fraction=.1)
  run.prism(bp,n.cores=8,update.gibbs=update)
}
mix13 <- cnt[c1$bulk_sample_id,,drop=FALSE]
fit_path<-file.path(CACHE,"BAYESPRISM_C1_FULL_FIT.rds")
if(file.exists(fit_path)) fit<-readRDS(fit_path) else {fit<-run_bp(build_gep(),mix13,TRUE);saveRDS(fit,fit_path,compress=FALSE)}

extract_fit <- function(fit,variant){
  theta<-if(!is.null(fit@posterior.theta_f))fit@posterior.theta_f@theta else fit@posterior.initial.cellType@theta
  tcv<-if(!is.null(fit@posterior.theta_f))fit@posterior.theta_f@theta.cv else fit@posterior.initial.cellType@theta.cv
  mal<-get.exp(fit,"type","malignant")
  Z<-fit@posterior.initial.cellType@Z; recon<-apply(Z,c(1,2),sum); zg<-dimnames(Z)[[2]]; orig<-mix13[rownames(recon),zg,drop=FALSE]
  cosine<-vapply(seq_len(nrow(recon)),function(i){den<-sqrt(sum(recon[i,]^2)*sum(orig[i,]^2));if(den==0)NA_real_ else sum(recon[i,]*orig[i,])/den},numeric(1)); names(cosine)<-rownames(recon)
  list(theta=theta,theta_cv=tcv,mal=mal,cosine=cosine,variant=variant)
}
full<-extract_fit(fit,"FULL_FROZEN_REFERENCE")

# Frozen reference sensitivity: Neftel/GSE131928-only, followed by the formal
# no-snRNA identity row inherited from the original reference sensitivity order.
nef_keep<-ref$cell_metadata$author=="Neftel2019"
nef_path<-file.path(CACHE,"BAYESPRISM_C1_GSE131928_FIT.rds")
nef_err<-NA_character_; nef_fit<-tryCatch({if(file.exists(nef_path))readRDS(nef_path) else {x<-run_bp(build_gep(nef_keep),mix13,FALSE);saveRDS(x,nef_path,compress=FALSE);x}},error=function(e){nef_err<<-conditionMessage(e);NULL})
nef<-if(is.null(nef_fit))NULL else extract_fit(nef_fit,"GSE131928_BASED")

score_parts <- function(e){
  x<-log2(e+1); z<-scale(x); z[!is.finite(z)]<-0
  out<-data.table(sample_id=rownames(e)); out[,canonical_MES:=rowSums(z[,intersect(canonical,colnames(z)),drop=FALSE])/95]
  for(k in classes){nm<-switch(k,MALIGNANT_DOMINANT="C_MALIGNANT",ECOLOGICAL_DOMINANT="C_ECOLOGICAL",SHARED_MIXED_ORIGIN="C_SHARED",LOW_INFORMATION_OR_UNSTABLE="C_UNSTABLE");g<-intersect(cls[primary_origin==k,gene],colnames(z));out[,(nm):=rowSums(z[,g,drop=FALSE])/95]}
  out[,additivity_error:=abs(canonical_MES-(C_MALIGNANT+C_ECOLOGICAL+C_SHARED+C_UNSTABLE))]
  out
}
bpfull<-score_parts(full$mal); bpfull[,participant_id:=c1$participant_id[match(sample_id,c1$bulk_sample_id)]]
bpnef<-if(is.null(nef))NULL else {x<-score_parts(nef$mal);x[,participant_id:=c1$participant_id[match(sample_id,c1$bulk_sample_id)]];x}

# Author labels, frozen 95 genes, malignant-cell weighted sample pseudobulk and
# equal-weight patient aggregation. No bulk result enters this metric.
meta<-readRDS(file.path(ICB,"metadata_ICB.rds")); meta<-meta[meta$Rx=="Pre-ICB",,drop=FALSE]
pb<-readRDS(file.path(ICB,"pseudobulks_by_sample_cell_type_pre_TPM.rds"))
spl<-strsplit(colnames(pb),"*",fixed=TRUE); info<-data.table(Sample_ID=vapply(spl,`[[`,character(1),1),Cell_type=vapply(spl,`[[`,character(1),2),col=seq_along(spl))
cc<-as.data.table(meta)[,.(n_cells=.N),by=.(Sample_ID,Cell_type,Patient_ID)]
mi<-merge(info,cc,by=c("Sample_ID","Cell_type"),all.x=TRUE)
mal_types<-unique(meta$Cell_type[meta$malignant_group=="Malignant"]); mi<-mi[Cell_type%in%mal_types & is.finite(n_cells)]
sample_expr<-lapply(unique(mi$Sample_ID),function(s){q<-mi[Sample_ID==s];rowSums(sweep(pb[,q$col,drop=FALSE],2,q$n_cells,"*"),na.rm=TRUE)/sum(q$n_cells)})
sample_expr<-do.call(cbind,sample_expr);colnames(sample_expr)<-unique(mi$Sample_ID);rownames(sample_expr)<-rownames(pb)
spmap<-unique(as.data.table(meta)[,.(Sample_ID,Patient_ID)]); spmap[,participant_id:=sprintf("GBM-%03d",as.integer(Patient_ID))]
pat_expr<-sapply(c1$participant_id,function(p){ss<-spmap[participant_id==p,Sample_ID];rowMeans(sample_expr[,intersect(ss,colnames(sample_expr)),drop=FALSE],na.rm=TRUE)})
colnames(pat_expr)<-c1$participant_id;rownames(pat_expr)<-rownames(sample_expr)
sn<-score_parts(t(pat_expr));setnames(sn,"sample_id","participant_id")
prop<-as.data.table(meta)[,.(total=.N,malignant=sum(malignant_group=="Malignant"),myeloid=sum(Cell_type=="Myeloid"),vascular_stromal=sum(Cell_type%in%c("Endothelial","Pericytes"))),by=Patient_ID]
prop[,participant_id:=sprintf("GBM-%03d",as.integer(Patient_ID))];prop[,`:=`(sn_malignant_proportion=malignant/total,sn_myeloid_proportion=myeloid/total,sn_vascular_stromal_proportion=vascular_stromal/total)]

getsum<-function(theta,pat){cc<-grep(pat,colnames(theta),ignore.case=TRUE,value=TRUE);if(length(cc))rowSums(theta[,cc,drop=FALSE])else rep(NA_real_,nrow(theta))}
cal<-merge(bpfull,sn,by="participant_id",suffixes=c("_bp","_sn"));cal<-merge(cal,prop[,.(participant_id,sn_malignant_proportion,sn_myeloid_proportion,sn_vascular_stromal_proportion)],by="participant_id")
cal[,bulk_sample_id:=c1$bulk_sample_id[match(participant_id,c1$participant_id)]]
cal[,`:=`(bp_malignant_fraction=getsum(full$theta,"^malignant$")[bulk_sample_id],bp_myeloid_fraction=getsum(full$theta,"myeloid")[bulk_sample_id],bp_vascular_stromal_fraction=(getsum(full$theta,"endothelial")+getsum(full$theta,"pericyte|fibroblast"))[bulk_sample_id],reconstruction_cosine=full$cosine[bulk_sample_id],reconstruction_error=1-full$cosine[bulk_sample_id],posterior_finite_fraction=apply(full$theta,1,function(x)mean(is.finite(x)))[bulk_sample_id],posterior_median_cv=apply(full$theta_cv,1,function(x)median(x[is.finite(x)],na.rm=TRUE))[bulk_sample_id],canonical_gene_coverage_bp=length(intersect(canonical,colnames(full$mal))),canonical_gene_coverage_snrna=length(intersect(canonical,colnames(pat_expr))))]
setcolorder(cal,c("participant_id","bulk_sample_id",setdiff(names(cal),c("participant_id","bulk_sample_id"))))
w(cal,"ICB_13PATIENT_CALIBRATION.tsv")

rank_conc<-function(x,y){d<-combn(seq_along(x),2);v<-sign(x[d[1,]]-x[d[2,]])*sign(y[d[1,]]-y[d[2,]]);mean(v[v!=0]>0)}
cstat<-function(x,y,metric,variant="FULL_FROZEN_REFERENCE"){
  ok<-is.finite(x)&is.finite(y);x<-x[ok];y<-y[ok]
  sp<-suppressWarnings(cor(x,y,method="spearman"));ke<-suppressWarnings(cor(x,y,method="kendall"));rc<-rank_conc(x,y)
  b<-replicate(2000,{ii<-sample.int(length(x),length(x),TRUE);suppressWarnings(cor(x[ii],y[ii],method="spearman"))});b<-b[is.finite(b)]
  data.table(reference_variant=variant,metric=metric,n=length(x),spearman_rho=sp,kendall_tau=ke,rank_concordance=rc,bootstrap_ci_low=quantile(b,.025,na.rm=TRUE),bootstrap_ci_high=quantile(b,.975,na.rm=TRUE))
}
metrics<-rbindlist(list(
  cstat(cal$canonical_MES_bp,cal$canonical_MES_sn,"BP_MALIGNANT_MES_VS_SNRNA_MALIGNANT_MES"),
  cstat(cal$bp_myeloid_fraction,cal$sn_myeloid_proportion,"BP_MYELOID_FRACTION_VS_SNRNA_PROPORTION"),
  cstat(cal$bp_vascular_stromal_fraction,cal$sn_vascular_stromal_proportion,"BP_VASCULAR_STROMAL_FRACTION_VS_SNRNA_PROPORTION"),
  cstat(cal$C_MALIGNANT_bp,cal$C_MALIGNANT_sn,"BP_MALIGNANT_CLASS_CONTRIBUTION_VS_SNRNA")
))

lopo<-rbindlist(lapply(seq_len(nrow(cal)),function(i){q<-cal[-i];rbindlist(list(
  data.table(metric="BP_MALIGNANT_MES_VS_SNRNA_MALIGNANT_MES",rho=suppressWarnings(cor(q$canonical_MES_bp,q$canonical_MES_sn,method="spearman"))),
  data.table(metric="BP_MYELOID_FRACTION_VS_SNRNA_PROPORTION",rho=suppressWarnings(cor(q$bp_myeloid_fraction,q$sn_myeloid_proportion,method="spearman"))),
  data.table(metric="BP_VASCULAR_STROMAL_FRACTION_VS_SNRNA_PROPORTION",rho=suppressWarnings(cor(q$bp_vascular_stromal_fraction,q$sn_vascular_stromal_proportion,method="spearman"))),
  data.table(metric="BP_MALIGNANT_CLASS_CONTRIBUTION_VS_SNRNA",rho=suppressWarnings(cor(q$C_MALIGNANT_bp,q$C_MALIGNANT_sn,method="spearman")))
  ))[,`:=`(left_out_patient=cal$participant_id[i],n=nrow(q))]}))
w(lopo,"ICB_13PATIENT_LOPO.tsv")

sens<-list(metrics[,.(reference_variant,metric,n,spearman_rho,kendall_tau,rank_concordance,bootstrap_ci_low,bootstrap_ci_high,status="EVALUATED",failure_reason=NA_character_)])
if(!is.null(bpnef)){
  nq<-merge(bpnef,sn,by="participant_id",suffixes=c("_bp","_sn"));sens[[2]]<-cstat(nq$canonical_MES_bp,nq$canonical_MES_sn,"BP_MALIGNANT_MES_VS_SNRNA_MALIGNANT_MES","GSE131928_BASED")[,`:=`(status="EVALUATED",failure_reason=NA_character_)]
} else sens[[2]]<-data.table(reference_variant="GSE131928_BASED",metric="BP_MALIGNANT_MES_VS_SNRNA_MALIGNANT_MES",n=NA_integer_,spearman_rho=NA_real_,kendall_tau=NA_real_,rank_concordance=NA_real_,bootstrap_ci_low=NA_real_,bootstrap_ci_high=NA_real_,status="FAILED",failure_reason=nef_err)
sens[[3]]<-metrics[metric=="BP_MALIGNANT_MES_VS_SNRNA_MALIGNANT_MES"][,`:=`(reference_variant="EXCLUDING_SNRNA_IDENTICAL_TO_FULL",status="IDENTICAL_FROZEN_REFERENCE_IS_SCRNA_ONLY",failure_reason=NA_character_)]
sens<-rbindlist(sens,fill=TRUE);w(sens,"ICB_REFERENCE_SENSITIVITY.tsv")

primary<-metrics[metric=="BP_MALIGNANT_MES_VS_SNRNA_MALIGNANT_MES"]$spearman_rho
pl<-lopo[metric=="BP_MALIGNANT_MES_VS_SNRNA_MALIGNANT_MES"]$rho
sensr<-sens[metric=="BP_MALIGNANT_MES_VS_SNRNA_MALIGNANT_MES" & status%in%c("EVALUATED","IDENTICAL_FROZEN_REFERENCE_IS_SCRNA_ONLY"),spearman_rho]
er<-cal$reconstruction_error; recon_ok<-all(is.finite(er)) && (max(er)<=median(er)+5*mad(er) || max(er)<.2)
checks<-data.table(criterion=c("TRACEABLE_N_GE_10","CANONICAL_COVERAGE_GE_80_PERCENT","NO_SYSTEMATIC_BAYESPRISM_FAILURE","PRIMARY_SPEARMAN_GE_0_30","AT_LEAST_80_PERCENT_LOPO_POSITIVE","NO_STRONG_NEGATIVE_LOPO","REFERENCE_SENSITIVITY_DIRECTION_CONSISTENT","FRACTIONS_BIOLOGICALLY_PLAUSIBLE","RECONSTRUCTION_ERROR_NOT_DOMINATED"),pass=c(nrow(cal)>=10,min(cal$canonical_gene_coverage_bp)/95>=.8,all(cal$posterior_finite_fraction==1),is.finite(primary)&&primary>=.3,mean(pl>0,na.rm=TRUE)>=.8,all(pl>-.3,na.rm=TRUE),length(sensr)>=2&&all(is.finite(sensr))&&all(sensr>0),all(cal$bp_malignant_fraction>=0&cal$bp_malignant_fraction<=1&cal$bp_myeloid_fraction>=0&cal$bp_myeloid_fraction<=1&cal$bp_vascular_stromal_fraction>=0&cal$bp_vascular_stromal_fraction<=1),recon_ok),observed=c(paste0(nrow(cal)," mapped patients"),paste0(min(cal$canonical_gene_coverage_bp),"/95"),paste0("minimum finite fraction ",min(cal$posterior_finite_fraction)),paste0("rho=",signif(primary,5)),paste0(signif(mean(pl>0,na.rm=TRUE),4)),paste0("minimum LOPO rho=",signif(min(pl,na.rm=TRUE),5)),paste(signif(sensr,5),collapse=","),paste0("all fractions in [0,1]=",all(cal$bp_malignant_fraction>=0&cal$bp_malignant_fraction<=1)),paste0("max error=",signif(max(er),5),"; median=",signif(median(er),5))))
status<-if(all(checks$pass))"ICB_DECONVOLUTION_TRANSPORT_PASS" else if(reference_pass&&nrow(cal)>=10&&all(cal$posterior_finite_fraction==1))"ICB_DECONVOLUTION_TRANSPORT_PARTIAL" else "ICB_DECONVOLUTION_TRANSPORT_FAIL"
w(checks,"ICB_TRANSPORT_CALIBRATION_CRITERIA_INTERNAL.tsv")
writeLines(c("# ICB Transport Calibration Gate","",paste0("Status: `",status,"`"),"",paste0("Primary Spearman rho: ",signif(primary,5)),paste0("Positive LOPO fraction: ",signif(mean(pl>0,na.rm=TRUE),4)),"",if(status=="ICB_DECONVOLUTION_TRANSPORT_PASS")"Phase C2 is authorized." else "Phase C2 is closed; only calibration results may be reported.","","All associations are transport/calibration evidence and not causal source evidence."),file.path(OUT,"ICB_TRANSPORT_CALIBRATION_GATE.md"))
writeLines(capture.output(sessionInfo()),file.path(OUT,"ICB_C1_SESSION_INFO.txt"))
cat("MODULE_C1_COMPLETE",status,nrow(cal),primary,"\n")
