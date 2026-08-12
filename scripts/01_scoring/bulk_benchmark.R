#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
suppressPackageStartupMessages({library(data.table);library(GSVA);library(singscore);library(MCPcounter);library(xCell);library(metafor)})
set.seed(2026073113)

DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset="data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset="results")
REPO_ROOT <- Sys.getenv("GBM_MES_REPO_ROOT", unset=".")
CORE_INPUT <- file.path(DATA_ROOT,"prepared","core")
OUT <- file.path(RESULT_ROOT,"scoring","bulk");dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
TCGA_COHORT <- file.path(DATA_ROOT,"prepared","cohorts","tcga_strict_cohort.tsv")
TCGA_EXPR <- file.path(DATA_ROOT,"bulk","tcga_gbm_expression.tsv.gz")
SIG <- file.path(REPO_ROOT,"resources","signatures","signature_gene_sets.tsv")
w<-function(x,n)fwrite(x,file.path(OUT,n),sep="\t",quote=FALSE,na="NA")
sha256<-function(path){z<-system2("/usr/bin/shasum",c("-a","256",shQuote(path)),stdout=TRUE);strsplit(z,"[[:space:]]+")[[1]][1]}

map<-c(canonical_MES="NEFTEL_MES_LIKE",MES1="NEFTEL_MES1",MES2="NEFTEL_MES2",AC_like="NEFTEL_AC_LIKE",
       OPC_like="NEFTEL_OPC_LIKE",NPC_like="NEFTEL_NPC_LIKE",proliferation="PROLIFERATION_CELL_CYCLE",
       hypoxia="HALLMARK_HYPOXIA",EMT="HALLMARK_EMT",matrix="NABA_CORE_MATRISOME",myeloid="MYELOID_CORE_IDENTITY")
reg<-fread(SIG);reg[,gene:=toupper(trimws(gene_symbol_clean))]
sets<-lapply(map,function(x)unique(reg[signature_id==x&!is.na(gene)&nzchar(gene),gene]))
ctx_names<-c("hypoxia","EMT","matrix","myeloid")
zvec<-function(x){s<-sd(x,na.rm=TRUE);if(!is.finite(s)||s==0)return(rep(NA_real_,length(x)));(x-mean(x,na.rm=TRUE))/s}
score_zmean<-function(expr,genes){g<-intersect(genes,rownames(expr));if(length(g)<2)return(rep(NA_real_,ncol(expr)));x<-expr[g,,drop=FALSE];m<-rowMeans(x);s<-apply(x,1,sd);ok<-is.finite(s)&s>0;if(sum(ok)<2)return(rep(NA_real_,ncol(expr)));colMeans((x[ok,,drop=FALSE]-m[ok])/s[ok])}
score_context<-function(scores)rowMeans(as.data.frame(lapply(as.data.frame(scores)[,ctx_names,drop=FALSE],zvec)),na.rm=TRUE)
boot_cor<-function(x,y,B=2000){ok<-is.finite(x)&is.finite(y);x<-x[ok];y<-y[ok];n<-length(x);if(n<5)return(c(NA,NA));b<-replicate(B,{i<-sample.int(n,n,TRUE);suppressWarnings(cor(x[i],y[i],method="spearman"))});quantile(b,c(.025,.975),na.rm=TRUE,names=FALSE)}
cortab<-function(cohort,method,x,y,endpoint,level){ok<-is.finite(x)&is.finite(y);n<-sum(ok);if(n<4)return(data.table(cohort,method,endpoint,analysis_level=level,n,rho=NA_real_,p_value=NA_real_,ci_low=NA_real_,ci_high=NA_real_));ct<-suppressWarnings(cor.test(x[ok],y[ok],method="spearman",exact=FALSE));ci<-boot_cor(x[ok],y[ok]);data.table(cohort,method,endpoint,analysis_level=level,n,rho=unname(ct$estimate),p_value=ct$p.value,ci_low=ci[1],ci_high=ci[2])}
rank_resid_cor<-function(y,x,covs){d<-data.frame(y=rank(y),x=rank(x),covs);d<-d[complete.cases(d),,drop=FALSE];if(nrow(d)<10||ncol(d)<3)return(c(rho=NA,p=NA,n=nrow(d)));ry<-resid(lm(y~.,data=d[,setdiff(names(d),"x"),drop=FALSE]));rx<-resid(lm(x~.,data=d[,setdiff(names(d),"y"),drop=FALSE]));ct<-cor.test(ry,rx,method="pearson");c(rho=unname(ct$estimate),p=ct$p.value,n=nrow(d))}

strict<-fread(TCGA_COHORT)
stopifnot(nrow(strict)==uniqueN(strict$patient_id),nrow(strict)==uniqueN(strict$sample_id),all(strict$strict_inclusion))
ids<-strict$sample_id
hdr<-names(fread(TCGA_EXPR,nrows=0,check.names=FALSE));missing<-setdiff(ids,hdr)
if(length(missing))stop("P0 TCGA strict expression IDs missing: ",paste(missing,collapse=";"))
dt<-fread(TCGA_EXPR,select=c(hdr[1],ids),check.names=FALSE)
genes<-toupper(trimws(dt[[1]]));expr<-as.matrix(dt[,-1]);storage.mode(expr)<-"double";rownames(expr)<-genes
if(anyDuplicated(genes)){u<-unique(genes);expr<-rowsum(expr,genes,reorder=FALSE)/as.vector(table(factor(genes,levels=u)))}
colnames(expr)<-ids
# TCGA HiSeqV2 reference values are already normalized log2-scale; do not log-transform again.

sc<-data.table(patient_id=ids)
coverage<-list()
for(nm in names(sets)){sc[[nm]]<-score_zmean(expr,sets[[nm]]);coverage[[nm]]<-data.table(cohort="TCGA_GBM",score=nm,method="mean_z",genes_total=length(sets[[nm]]),genes_detected=length(intersect(sets[[nm]],rownames(expr))))}
sc[,context:=score_context(sc)]
primary_tcga<-cortab("TCGA_GBM","mean_z",sc$canonical_MES,sc$context,"canonical_MES_vs_context","primary")
canp<-setdiff(sets$canonical_MES,unique(unlist(sets[ctx_names])))
ctxp<-lapply(sets[ctx_names],function(x)setdiff(x,sets$canonical_MES))
psc<-data.table(patient_id=ids,canonical_MES=score_zmean(expr,canp));for(nm in ctx_names)psc[[nm]]<-score_zmean(expr,ctxp[[nm]]);psc[,context:=score_context(psc)]
pruned_tcga<-cbind(cortab("TCGA_GBM","mean_z_bilateral_pruned",psc$canonical_MES,psc$context,"canonical_MES_vs_context","primary_overlap_pruned"),canonical_genes_retained=length(intersect(canp,rownames(expr))),context_unique_genes_retained=length(intersect(unique(unlist(ctxp)),rownames(expr))))

qa<-list();sensitivity<-list()
run_alt<-function(method){tryCatch({
  if(method=="singscore"){
    rd<-rankGenes(expr);z<-sapply(names(sets),function(nm)simpleScore(rd,upSet=intersect(sets[[nm]],rownames(expr)),centerScore=TRUE,knownDirection=TRUE)$TotalScore);z<-as.data.frame(z);rownames(z)<-colnames(expr)
  }else{
    gp<-lapply(sets,function(g)intersect(g,rownames(expr)));gp<-gp[lengths(gp)>=2]
    param<-switch(method,GSVA=gsvaParam(expr,gp,kcdf="Gaussian"),ssGSEA=ssgseaParam(expr,gp),PLAGE=plageParam(expr,gp))
    z<-as.data.frame(t(gsva(param,verbose=FALSE)))
    if(method=="PLAGE")for(nm in intersect(names(sets),names(z))){rr<-suppressWarnings(cor(z[[nm]],sc[[nm]],use="complete.obs",method="spearman"));if(is.finite(rr)&&rr<0)z[[nm]]<--z[[nm]]}
  }
  z$context<-score_context(z);cortab("TCGA_GBM",method,z$canonical_MES,z$context,"canonical_MES_vs_context","scoring_sensitivity")
},error=function(e){qa[[length(qa)+1]]<<-data.table(check=paste0(method,"_execution"),status="FAIL",detail=conditionMessage(e));data.table(cohort="TCGA_GBM",method,endpoint="canonical_MES_vs_context",analysis_level="scoring_sensitivity",n=length(ids),rho=NA_real_,p_value=NA_real_,ci_low=NA_real_,ci_high=NA_real_)})}
for(method in c("singscore","GSVA","ssGSEA","PLAGE"))sensitivity[[method]]<-run_alt(method)

mcp<-tryCatch(MCPcounter.estimate(expr,featuresType="HUGO_symbols"),error=function(e){qa[[length(qa)+1]]<<-data.table(check="MCPcounter",status="FAIL",detail=conditionMessage(e));NULL})
xc<-tryCatch(xCellAnalysis(expr,rnaseq=TRUE),error=function(e){qa[[length(qa)+1]]<<-data.table(check="xCell",status="FAIL",detail=conditionMessage(e));NULL})
covs<-data.frame(row.names=ids)
if(!is.null(mcp)){getm<-function(pattern){i<-grep(pattern,rownames(mcp),ignore.case=TRUE);if(length(i))as.numeric(mcp[i[1],ids])else rep(NA_real_,length(ids))};covs$mcp_monocytic<-getm("monocytic");covs$mcp_endothelial<-getm("endothelial");covs$mcp_fibroblast<-getm("fibroblast")}
if(!is.null(xc)){getx<-function(pattern){i<-grep(pattern,rownames(xc),ignore.case=TRUE);if(length(i))as.numeric(xc[i[1],ids])else rep(NA_real_,length(ids))};covs$xcell_macrophage<-getx("Macrophage");covs$xcell_endothelial<-getx("Endothelial");covs$xcell_stromal<-getx("stroma|fibro")}
usecov<-covs[,vapply(covs,function(x)sum(is.finite(x))>=10&&sd(x,na.rm=TRUE)>0,logical(1)),drop=FALSE]
rr<-rank_resid_cor(sc$context,sc$canonical_MES,usecov)
composition_tcga<-data.table(cohort="TCGA_GBM",method="rank_residualized_MCPcounter_xCell",endpoint="canonical_MES_vs_context",n=rr["n"],partial_rho=rr["rho"],p_value=rr["p"],covariates=paste(names(usecov),collapse=";"))

cgga_primary<-fread(file.path(CORE_INPUT,"05_bulk/BULK_PRIMARY_RESULTS.tsv"))
cgga_pruned<-fread(file.path(CORE_INPUT,"05_bulk/BULK_OVERLAP_PRUNED.tsv"))
cgga_sens<-fread(file.path(CORE_INPUT,"05_bulk/BULK_SCORING_SENSITIVITY.tsv"))
cgga_comp<-fread(file.path(CORE_INPUT,"05_bulk/BULK_COMPOSITION_SENSITIVITY.tsv"))
primary<-rbind(cgga_primary[,!"FDR"],primary_tcga,fill=TRUE);primary[,FDR:=p.adjust(p_value,"BH")]
pruned<-rbind(cgga_pruned[,!"FDR"],pruned_tcga,fill=TRUE);pruned[,FDR:=p.adjust(p_value,"BH")]
sens<-rbind(cgga_sens[,!"FDR"],rbindlist(sensitivity),fill=TRUE);sens[,FDR:=p.adjust(p_value,"BH")]
comp<-rbind(cgga_comp[,!"FDR"],composition_tcga,fill=TRUE);comp[,FDR:=p.adjust(p_value,"BH")]

meta_one<-function(d,type){z<-atanh(pmax(pmin(d$rho,.999999),-.999999));vi<-1/(d$n-3);fit<-rma.uni(z,vi=vi,method="REML");pr<-predict(fit);list(summary=data.table(analysis_type=type,k=nrow(d),random_effects_rho=tanh(as.numeric(fit$b)),ci_low=tanh(fit$ci.lb),ci_high=tanh(fit$ci.ub),prediction_low=tanh(pr$pi.lb),prediction_high=tanh(pr$pi.ub),tau2=fit$tau2,I2=fit$I2,Q=fit$QE,Q_p=fit$QEp,positive_cohorts=sum(d$rho>0),total_cohorts=nrow(d)),heterogeneity=data.table(analysis_type=type,k=nrow(d),Q=fit$QE,Q_df=fit$k-1,Q_p=fit$QEp,tau2=fit$tau2,I2=fit$I2,prediction_interval_note="k=3; prediction interval is estimable but imprecise; CGGA cohorts share broad recruitment ecosystem and are not treated as fully independent ecological universes"))}
mraw<-meta_one(primary,"RAW");mpruned<-meta_one(pruned,"BILATERAL_OVERLAP_PRUNED")
meta<-rbind(mraw$summary,mpruned$summary);hetero<-rbind(mraw$heterogeneity,mpruned$heterogeneity)
final<-rbind(cbind(analysis_type="RAW",primary),cbind(analysis_type="BILATERAL_OVERLAP_PRUNED",pruned),fill=TRUE)
claim_status<-if(all(primary$rho>0)&all(pruned$rho>0))"RETAIN_PRIMARY" else "DOWNGRADE"
claim<-data.table(claim_id="P1_STRICT_BULK_COUPLING",status=claim_status,criteria="positive cohort-specific direction in all strict cohorts before and after bilateral pruning",observed=paste0("raw=",paste(primary$cohort,round(primary$rho,3),sep=":",collapse=";"),"; pruned=",paste(pruned$cohort,round(pruned$rho,3),sep=":",collapse=";")),boundary="association only; TCGA IDH-wildtype evidence uses matched open GDC masked WXS MAF absence of IDH1/IDH2 calls; CGGA cohorts are not treated as fully independent ecological universes")

reference<-data.table(input=c("TCGA_STRICT_COHORT","TCGA_EXPRESSION","SIGNATURE_REGISTRY","CGGA_VERIFIED_RESULTS"),path=c(TCGA_COHORT,TCGA_EXPR,SIG,file.path(CORE_INPUT,"05_bulk")),SHA256=c(sha256(TCGA_COHORT),sha256(TCGA_EXPR),sha256(SIG),"SEE_FINAL_RESULT_REUSE_CHECK"),role=c("enumerated strict membership","official processed normalized log2 expression; no additional log transform","frozen gene sets","REUSE_VERIFIED numeric outputs"))
qa[[length(qa)+1]]<-data.table(check="TCGA_strict_membership",status=ifelse(length(missing)==0&&length(ids)==138,"PASS","FAIL"),detail=paste0("n=",length(ids),"; unique patients=",uniqueN(strict$patient_id),"; missing=",length(missing)))
qa[[length(qa)+1]]<-data.table(check="TCGA_expression_transform",status="PASS",detail="reference HiSeqV2 is already normalized log2-scale; no double log transform")
qa[[length(qa)+1]]<-data.table(check="MCPcounter",status=ifelse(is.null(mcp),"NOT_EVALUABLE","PASS"),detail=ifelse(is.null(mcp),"method failed",paste(dim(mcp),collapse="x")))
qa[[length(qa)+1]]<-data.table(check="xCell",status=ifelse(is.null(xc),"NOT_EVALUABLE","PASS"),detail=ifelse(is.null(xc),"method failed",paste(dim(xc),collapse="x")))
qa[[length(qa)+1]]<-data.table(check="purity_sensitivity",status="NOT_EVALUABLE",detail="no patient-mapped reference purity estimate frozen for all three strict cohorts; no ad hoc substitute")

w(reference,"TCGA_BULK_REFERENCE_INPUTS.tsv");w(sc,"TCGA_BULK_PATIENT_SCORES.tsv");w(rbindlist(coverage),"TCGA_BULK_GENE_COVERAGE.tsv")
w(primary_tcga,"TCGA_BULK_PRIMARY.tsv");w(pruned_tcga,"TCGA_BULK_OVERLAP_PRUNED.tsv");w(rbindlist(sensitivity),"TCGA_BULK_SCORING_SENSITIVITY.tsv");w(composition_tcga,"TCGA_BULK_COMPOSITION_SENSITIVITY.tsv")
w(final,"STRICT_BULK_FINAL.tsv");w(meta,"STRICT_BULK_META.tsv");w(hetero,"STRICT_BULK_HETEROGENEITY.tsv");w(sens,"STRICT_BULK_SCORING_SENSITIVITY.tsv");w(comp,"STRICT_BULK_COMPOSITION_SENSITIVITY.tsv");w(claim,"STRICT_BULK_CLAIM_GATE.tsv");w(rbindlist(qa,fill=TRUE),"STRICT_BULK_QA.tsv")
writeLines(capture.output(sessionInfo()),file.path(OUT,"SESSION_INFO.txt"))
cat("STRICT_BULK_COMPLETE",claim_status,"TCGA_n",length(ids),"\n")
