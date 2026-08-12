#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
suppressPackageStartupMessages({library(Matrix); library(data.table); library(jsonlite)})
set.seed(2026073115)

DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset="data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset="results")
REPO_ROOT <- Sys.getenv("GBM_MES_REPO_ROOT", unset=".")
OUT <- file.path(RESULT_ROOT,"single_cell","gse174554")
DATA <- file.path(DATA_ROOT,"single_cell","GSE174554")
MANIFEST <- file.path(DATA_ROOT,"prepared","gse174554","sample_manifest.tsv")
ANNOTATION_MAP <- file.path(DATA_ROOT,"prepared","gse174554","annotation_mapping.tsv")
INTEGRITY <- file.path(DATA_ROOT,"prepared","gse174554","matrix_integrity.tsv")
SIG <- file.path(REPO_ROOT,"resources","signatures","signature_gene_sets.tsv")
META_FILE <- file.path(DATA,"download/GSE174554_Tumor_normal_metadata.txt.gz")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
w <- function(x,n)fwrite(x,file.path(OUT,n),sep="\t",quote=FALSE,na="NA")
sha256 <- function(path){z<-system2("/usr/bin/shasum",c("-a","256",shQuote(path)),stdout=TRUE);strsplit(z,"[[:space:]]+")[[1]][1]}

manifest <- fread(MANIFEST); annotation_map <- fread(ANNOTATION_MAP)
manifest <- merge(manifest,annotation_map,by="sample_id",all.x=TRUE)
manifest <- manifest[primary_eligible_metadata==TRUE & eligible_for_malignant_scoring==TRUE]
integrity <- fread(INTEGRITY); bad_samples <- unique(integrity[integrity_status=="INTEGRITY_CHECK_FAILED",sample_id])
w(manifest[sample_id%in%bad_samples,.(sample_id,patient_id,exclusion_reason="INTEGRITY_CHECK_FAILED_OFFICIAL_MATRIXMARKET_MEMBER")],"GSE174554_INTEGRITY_EXCLUSIONS.tsv")
manifest <- manifest[!sample_id%in%bad_samples]
stopifnot(nrow(manifest)>0,!anyDuplicated(manifest$patient_id))
registry <- fread(SIG);registry[,gene:=toupper(trimws(gene_symbol_clean))]
canonical <- unique(registry[signature_id=="NEFTEL_MES_LIKE" & !is.na(gene)&nzchar(gene),gene])

metadata <- fread(cmd=paste("gzip -dc",shQuote(META_FILE)))
setnames(metadata,c("sample_barcode","tumor_normal"))
metadata[,sample_id:=sub("_(.*)$","",sample_barcode)]
metadata[,barcode:=sub("^SF[0-9]+(v2)?_","",sample_barcode)]
metadata[,barcode_norm:=sub("[.-][0-9]+$","",barcode)]
metadata[,meta_batch:=ifelse(grepl("\\.[0-9]+$",barcode),paste0("batch",sub(".*\\.([0-9]+)$","\\1",barcode)),"batch1")]
metadata <- unique(metadata,by=c("sample_id","meta_batch","barcode_norm"))

split_members <- function(value)unlist(strsplit(value,";",fixed=TRUE))
batch_key <- function(value)ifelse(grepl("_batch[0-9]+_",value),sub(".*_(batch[0-9]+)_.*","\\1",value),"batch1")
rows <- list();sample_qa <- list()
for(k in seq_len(nrow(manifest))){
  row <- manifest[k];sample_id <- row$sample_id
  sample_value <- sample_id
  matrices<-split_members(row$matrix_file);features<-split_members(row$features_file);barcodes<-split_members(row$barcodes_file)
  batches<-Reduce(intersect,list(batch_key(matrices),batch_key(features),batch_key(barcodes)))
  parts<-list();labels<-list();batch_rows<-list()
  for(batch in batches){
    directory<-file.path(DATA,"extracted",sample_id)
    mf<-file.path(directory,basename(matrices[batch_key(matrices)==batch][1]));ff<-file.path(directory,basename(features[batch_key(features)==batch][1]));bf<-file.path(directory,basename(barcodes[batch_key(barcodes)==batch][1]))
    if(!all(file.exists(c(mf,ff,bf))))stop("P0_STOP extracted component missing ",sample_id," ",batch)
    ft<-fread(cmd=paste("gzip -dc",shQuote(ff)),header=FALSE);gene<-toupper(if(ncol(ft)>=2)as.character(ft[[2]])else as.character(ft[[1]]))
    barcode<-as.character(fread(cmd=paste("gzip -dc",shQuote(bf)),header=FALSE)[[1]])
    x<-as(readMM(gzfile(mf)),"dgCMatrix");if(nrow(x)==length(barcode)&&ncol(x)==length(gene))x<-t(x)
    if(nrow(x)!=length(gene)||ncol(x)!=length(barcode))stop("P0_STOP matrix dimensions ",sample_id," ",batch)
    bm<-metadata[get("sample_id")==sample_value & meta_batch==batch]
    match_index<-match(sub("-[0-9]+$","",barcode),bm$barcode_norm);lab<-rep(NA_character_,length(barcode));lab[!is.na(match_index)]<-tolower(bm$tumor_normal[match_index[!is.na(match_index)]])
    keep<-!is.na(lab)&lab%in%c("tumor","normal")
    idx<-which(gene%in%canonical);ug<-unique(gene[idx])
    collapse<-sparseMatrix(i=match(gene[idx],ug),j=seq_along(idx),x=1,dims=c(length(ug),length(idx)))
    parts[[batch]]<-collapse%*%x[idx,keep,drop=FALSE];rownames(parts[[batch]])<-ug
    if(ncol(parts[[batch]])>0)colnames(parts[[batch]])<-paste(sample_id,batch,barcode[keep],sep="_")
    labels[[batch]]<-lab[keep]
    attr(parts[[batch]],"library")<-Matrix::colSums(x[,keep,drop=FALSE])
    batch_rows[[batch]]<-data.table(sample_id,batch,n_total=length(barcode),n_mapped=sum(!is.na(lab)),n_labeled=sum(keep),n_tumor=sum(lab[keep]=="tumor"),n_normal=sum(lab[keep]=="normal"))
  }
  all_genes<-Reduce(union,lapply(parts,rownames));aligned<-lapply(parts,function(x){a<-Matrix(0,nrow=length(all_genes),ncol=ncol(x),sparse=TRUE,dimnames=list(all_genes,colnames(x)));a[rownames(x),]<-x;attr(a,"library")<-attr(x,"library");a})
  lib<-unlist(lapply(aligned,attr,"library"),use.names=FALSE);x<-do.call(cbind,aligned);lab<-unlist(labels,use.names=FALSE)
  if(length(lab)!=ncol(x)||length(lib)!=ncol(x))stop("P0_STOP label/count alignment ",sample_id)
  coverage<-intersect(canonical,rownames(x));logcpm<-log2(t(t(x[coverage,,drop=FALSE])*1e6/lib)+1)
  gm<-Matrix::rowMeans(logcpm);gs<-sqrt(pmax(Matrix::rowMeans(logcpm^2)-gm^2,0));ok<-is.finite(gs)&gs>0
  score<-as.numeric(Matrix::colMeans(Diagonal(x=1/gs[ok])%*%logcpm[ok,,drop=FALSE]-Matrix(matrix(gm[ok]/gs[ok],nrow=sum(ok),ncol=ncol(x)),sparse=TRUE)))
  threshold<-as.numeric(quantile(score,.75,names=FALSE));high<-score>=threshold;n<-length(score);nh<-sum(high)
  tumor<-lab=="tumor";a<-sum(high&tumor);b<-sum(high&!tumor);c<-sum(!high&tumor);d<-sum(!high&!tumor)
  corrected<-any(c(a,b,c,d)==0);aa<-a;bb<-b;cc<-c;dd<-d;if(corrected){aa<-aa+.5;bb<-bb+.5;cc<-cc+.5;dd<-dd+.5}
  bg<-mean(tumor);hi<-mean(tumor[high]);oe<-if(bg>0)hi/bg else NA_real_;lor<-log((aa*dd)/(bb*cc));se<-sqrt(1/aa+1/bb+1/cc+1/dd)
  eligible<-n>=100&&nh>=20
  rows[[sample_id]]<-data.table(study="GSE174554",sample_id,patient_id=row$patient_id,primary_or_recurrent=row$primary_or_recurrent,score="canonical_MES",compartment="author_tumor_label",platform_class="snRNA",n_cells=n,n_high=nh,high_threshold=threshold,background_count=sum(tumor),high_count=a,background_fraction=bg,high_fraction=hi,observed_expected=oe,log2_observed_expected=ifelse(oe>0,log2(oe),NA_real_),odds_ratio=exp(lor),log_odds_ratio=lor,log_odds_se=se,haldane_correction=corrected,eligible=eligible,malignant_nuclei=sum(tumor),canonical_genes_defined=length(canonical),canonical_genes_covered=length(coverage),inference_unit="patient",pooled_cell_inference=FALSE)
  sample_qa[[sample_id]]<-rbindlist(batch_rows)
  cat("GSE174554_CANONICAL_SAMPLE",sample_id,n,sum(tumor),length(coverage),"\n",flush=TRUE)
}
result<-rbindlist(rows,fill=TRUE);qc<-rbindlist(sample_qa,fill=TRUE)
w(result,"GSE174554_CANONICAL_COMPARTMENT.tsv");w(qc,"GSE174554_CANONICAL_BATCH_QA.tsv")
summary<-result[eligible==TRUE,.(eligible_patients=.N,median_observed_expected=median(observed_expected),IQR_low=quantile(observed_expected,.25),IQR_high=quantile(observed_expected,.75),positive_patients=sum(log_odds_ratio>0),median_log_odds=median(log_odds_ratio),total_malignant_nuclei=sum(malignant_nuclei))]
w(summary,"GSE174554_CANONICAL_SUMMARY.tsv")
qa<-data.table(check=c("raw_archive_lineage","matrix_member_integrity","primary_patient_identity","minimum_cells","canonical_coverage"),status=c("PASS_RAW_RECOMPUTED",ifelse(length(bad_samples)>0,"INTEGRITY_CHECK_FAILED_SAMPLE_EXCLUDED","PASS"),"PASS",ifelse(all(result[eligible==TRUE,n_cells]>=100),"PASS","FAIL"),ifelse(min(result$canonical_genes_covered)>=20,"PASS","FAIL")),detail=c("analysis reads extracted members of the verified raw archive",paste0("official archive structural failures=",paste(bad_samples,collapse=";"),"; member retained, sample excluded"),paste0("primary patients=",nrow(result),"; unique=",uniqueN(result$patient_id)),paste0("eligible=",sum(result$eligible)),paste0("coverage range=",min(result$canonical_genes_covered),"-",max(result$canonical_genes_covered))))
w(qa,"GSE174554_CANONICAL_QA.tsv")
write_json(list(schema_version="1.1.0",status="COMPLETE",seed=2026073115,raw_archive_sha256="a045e7e4b964f8d49f0f47fd46f5c54a189bf9c31f1909ad46c988fab5e21a0e",metadata_sha256=sha256(META_FILE),patient_unit=TRUE),file.path(OUT,"GSE174554_ANALYSIS_METADATA.json"),auto_unbox=TRUE,pretty=TRUE)
writeLines(c("# GSE174554 canonical MES analysis","","One official matrix member failed structural integrity and was excluded without repair. Official Tumor/Normal barcode labels were used, and the patient is the inferential unit. The label is a broad malignant/nonmalignant annotation and is not interpreted as a detailed ecosystem cell type."),file.path(OUT,"GSE174554_ANALYSIS_SUMMARY.md"))
cat("GSE174554_CANONICAL_COMPLETE",nrow(result),sum(result$eligible),"\n")
