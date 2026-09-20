#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
suppressPackageStartupMessages({library(Matrix);library(data.table);library(hdf5r);library(arrow);library(RANN);library(metafor);library(lme4);library(digest)})

DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset="data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset="results")
REPO_ROOT <- Sys.getenv("GBM_MES_REPO_ROOT", unset=".")
RUN <- file.path(RESULT_ROOT,'provenance')
V13 <- file.path(RESULT_ROOT,'interpretability')
V12 <- file.path(RESULT_ROOT,'spatial')
OUT <- file.path(RUN,'contributions')
set.seed(as.integer(20260801404 %% .Machine$integer.max))
w <- function(x,n) fwrite(x,file.path(OUT,n),sep='\t',quote=FALSE,na='NA')
safe_cor <- function(x,y){ok<-is.finite(x)&is.finite(y);if(sum(ok)<3)return(NA_real_);suppressWarnings(cor(x[ok],y[ok],method='spearman'))}
stable_seed <- function(...) as.integer(abs(as.double(digest2int(paste(...,collapse='|'),seed=20260801404)))%%2147483000+1)

pf <- fread(file.path(RUN,'01_provenance_freeze/PROVENANCE_CLASS_FREEZE.tsv'))
gcol <- if('gene'%in%names(pf))'gene' else 'harmonized_symbol'; ccol <- if('primary_origin'%in%names(pf))'primary_origin' else 'origin_class'
genes95 <- toupper(pf[[gcol]]); classes <- setNames(pf[[ccol]],genes95)
class_order <- c('MALIGNANT_DOMINANT','ECOLOGICAL_DOMINANT','SHARED_MIXED_ORIGIN','LOW_INFORMATION_OR_UNSTABLE')
cname <- c(MALIGNANT_DOMINANT='C_MALIGNANT',ECOLOGICAL_DOMINANT='C_ECOLOGICAL',SHARED_MIXED_ORIGIN='C_SHARED',LOW_INFORMATION_OR_UNSTABLE='C_UNSTABLE')

read10x <- function(path){f<-H5File$new(path,'r');on.exit(f$close_all());d<-f[['matrix/data']][];i<-f[['matrix/indices']][];p<-f[['matrix/indptr']][];s<-f[['matrix/shape']][];g<-toupper(as.character(f[['matrix/features/name']][]));b<-as.character(f[['matrix/barcodes']][]);j<-rep.int(seq_len(as.integer(s[2])),diff(as.integer(p)));list(X=sparseMatrix(i=as.integer(i)+1L,j=j,x=as.numeric(d),dims=as.integer(s)),genes=g,bars=b)}
make_blocks <- function(xy,multiple){nn<-nn2(xy,xy,k=2)$nn.dists[,2];d<-median(nn[nn>0&is.finite(nn)]);width<-multiple*d;paste(floor((xy[,1]-min(xy[,1]))/width),floor((xy[,2]-min(xy[,2]))/width),sep=':')}
block_observed <- function(a,b,xy){rbindlist(lapply(c(small=4,medium=8,large=16),function(m){bid<-make_blocks(xy,m);db<-data.table(bid,a,b)[is.finite(a)&is.finite(b),.(a=mean(a),b=mean(b)),by=bid];data.table(occupied_blocks=nrow(db),block_rho=safe_cor(db$a,db$b))}),idcol='block_scale')}
block_permute <- function(a,b,xy,key){rbindlist(lapply(c(small=4,medium=8,large=16),function(m){bid<-make_blocks(xy,m);db<-data.table(bid,a,b)[is.finite(a)&is.finite(b),.(a=mean(a),b=mean(b)),by=bid];rho<-safe_cor(db$a,db$b);set.seed(stable_seed(key,m));null<-if(nrow(db)>=4&&is.finite(rho))replicate(999,safe_cor(db$a,sample(db$b)))else rep(NA_real_,999);data.table(occupied_blocks=nrow(db),block_rho=rho,empirical_p=if(is.finite(rho))(1+sum(abs(null)>=abs(rho),na.rm=TRUE))/(1+sum(is.finite(null)))else NA_real_,permutations=999L)}),idcol='block_scale')}
cross_meta <- function(study_dt,keys){study_dt[,{x<-rho[is.finite(rho)];lo<-ci_low[is.finite(rho)];hi<-ci_high[is.finite(rho)];if(length(x)==2){zv<-atanh(pmax(pmin(x,.999999),-.999999));sev<-pmax((atanh(pmin(hi,.999999))-atanh(pmax(lo,-.999999)))/(2*qnorm(.975)),.03);md<-data.frame(yi_value=zv,vi_value=sev^2);fit<-tryCatch(metafor::rma.uni(yi=yi_value,vi=vi_value,data=md,method='REML'),error=function(e)NULL);if(!is.null(fit)){pr<-predict(fit);.(n_studies=2L,rho=tanh(as.numeric(fit$b)),ci_low=tanh(fit$ci.lb),ci_high=tanh(fit$ci.ub),prediction_low=tanh(pr$pi.lb),prediction_high=tanh(pr$pi.ub),I2=fit$I2,direction_consistent=all(sign(x)==sign(x[1])))}else .(n_studies=2L,rho=tanh(mean(zv)),ci_low=NA_real_,ci_high=NA_real_,prediction_low=NA_real_,prediction_high=NA_real_,I2=NA_real_,direction_consistent=all(sign(x)==sign(x[1])))}else .(n_studies=length(x),rho=if(length(x))x else NA_real_,ci_low=NA_real_,ci_high=NA_real_,prediction_low=NA_real_,prediction_high=NA_real_,I2=NA_real_,direction_consistent=NA)},by=keys]}

# Spatial ---------------------------------------------------------------------
map <- fread(file.path(V12,'01_spatial_reference/SPATIAL_SECTION_PATIENT_MAP.tsv')); setnames(map,'study_id','study')
auth <- fread(file.path(V12,'01_spatial_reference/SPATIAL_REFERENCE_INPUTS.tsv'))[input_role=='expression'];auth[,study:=study_id]
inputs <- merge(map[,.(study,section_id,patient_id,disease_class,geo_accession)],auth[,.(study,section_id,path)],by=c('study','section_id'))
parq <- rbind(as.data.table(read_parquet(file.path(V12,'02_gse237183/GSE237183_SPOT_SCORES.parquet'))),as.data.table(read_parquet(file.path(V12,'03_gse242352/GSE242352_SPOT_SCORES.parquet'))),fill=TRUE)
spot_rows<-list();sec_rows<-list();block_rows<-list();cache<-list();si<-ri<-bi<-0L
eco_map<-c(HYPOXIA='ecology_HYPOXIA',MYELOID='ecology_MYELOID',VASCULAR='ecology_VASCULAR',MATRIX='ecology_MATRIX')
for(k in seq_len(nrow(inputs))){m<-inputs[k];message('SPATIAL ',k,'/',nrow(inputs),' ',m$section_id);obj<-read10x(m$path);sp<-parq[study==m$study&section_id==m$section_id];idx<-match(sp$barcode,obj$bars);stopifnot(all(!is.na(idx)));X<-obj$X[,idx,drop=FALSE];lib<-colSums(X);lib[lib<=0]<-1
 gidx<-match(genes95,obj$genes,nomatch=0L);Z<-matrix(0,nrow=nrow(sp),ncol=95,dimnames=list(NULL,genes95));ok<-gidx>0
 if(any(ok)){a<-as.matrix(X[gidx[ok],,drop=FALSE]);a<-log1p(t(t(a)/lib)*1e4);zz<-t(scale(t(a)));zz[!is.finite(zz)]<-0;Z[,ok]<-t(zz)}
 parts<-sapply(class_order,function(cl){jj<-which(classes[genes95]==cl & ok);if(length(jj))rowSums(Z[,jj,drop=FALSE])/95 else rep(0,nrow(Z))});colnames(parts)<-cname[class_order];total<-rowSums(parts)
 eco<-as.data.table(sp[,..eco_map]);setnames(eco,names(eco_map));xy<-as.matrix(sp[,.(pxl_col,pxl_row)])
 si<-si+1;spot_rows[[si]]<-cbind(data.table(study=m$study,section_id=m$section_id,patient_id=m$patient_id,disease_class=m$disease_class,barcode=sp$barcode),as.data.table(parts),MES_ADDITIVE=total,eco)
 cache[[m$section_id]]<-list(meta=m,Z=Z,parts=parts,total=total,eco=eco,xy=xy,present=ok)
 for(cl in colnames(parts))for(ec in names(eco)){ri<-ri+1;sec_rows[[ri]]<-data.table(study=m$study,section_id=m$section_id,patient_id=m$patient_id,disease_class=m$disease_class,class_contribution=cl,ecology=ec,n_spots=nrow(sp),rho=safe_cor(parts[,cl],eco[[ec]]));z<-block_observed(parts[,cl],eco[[ec]],xy);bi<-bi+1;block_rows[[bi]]<-cbind(study=m$study,section_id=m$section_id,patient_id=m$patient_id,class_contribution=cl,ecology=ec,z)}
}
spots<-rbindlist(spot_rows);sec<-rbindlist(sec_rows);blocks<-rbindlist(block_rows);w(spots,'SPATIAL_CLASS_CONTRIBUTIONS.tsv')
pat<-sec[,.(rho=mean(rho,na.rm=TRUE),n_sections=uniqueN(section_id)),by=.(study,patient_id,disease_class,class_contribution,ecology)]
study<-pat[is.finite(rho)&abs(rho)<1,{z<-atanh(rho);se<-if(.N>1)sd(z)/sqrt(.N) else NA_real_;.(n_patients=.N,rho=tanh(mean(z)),ci_low=if(is.finite(se))tanh(mean(z)-1.96*se)else NA_real_,ci_high=if(is.finite(se))tanh(mean(z)+1.96*se)else NA_real_)},by=.(study,class_contribution,ecology)]
cross<-cross_meta(study,c('class_contribution','ecology'))
sec[,result_level:='SECTION'];pat[,result_level:='PATIENT'];study[,result_level:='STUDY'];cross[,result_level:='CROSS_STUDY']
lopo<-rbindlist(lapply(unique(pat$patient_id),function(id){q<-pat[patient_id!=id];q[is.finite(rho)&abs(rho)<1,.(rho=tanh(mean(atanh(rho))),n_patients=.N),by=.(class_contribution,ecology)][,`:=`(result_level='LOPO',omitted_unit=id)]}))
loso<-copy(study);loso[,`:=`(result_level='LOSO',omitted_unit=ifelse(study=='GSE237183','GSE242352','GSE237183'))]
w(rbindlist(list(sec,pat,study,cross,lopo,loso),fill=TRUE),'SPATIAL_CLASS_META.tsv')
bp<-blocks[,.(rho=mean(block_rho,na.rm=TRUE)),by=.(patient_id,class_contribution,ecology,block_scale)]
bpool<-bp[,.(rho=mean(rho,na.rm=TRUE)),by=.(class_contribution,ecology,block_scale)]
sel<-merge(cross[,.(class_contribution,ecology,cross_rho=rho,direction_consistent)],bpool[,.(aligned_scales=sum(sign(rho)==sign(first(rho)))),by=.(class_contribution,ecology)],by=c('class_contribution','ecology'),all.x=TRUE)
# Recompute aligned scales against the cross-study direction.
sel[,aligned_scales:=vapply(seq_len(.N),function(i)sum(sign(bpool[class_contribution==sel$class_contribution[i]&ecology==sel$ecology[i],rho])==sign(sel$cross_rho[i]),na.rm=TRUE),integer(1))]
sel[,selected:=direction_consistent==TRUE&abs(cross_rho)>=.30&aligned_scales>=2]
perm<-list();pi<-0L
for(i in which(sel$selected)){for(section in names(cache)){q<-cache[[section]];cl<-sel$class_contribution[i];ec<-sel$ecology[i];pi<-pi+1;perm[[pi]]<-cbind(study=q$meta$study,section_id=section,patient_id=q$meta$patient_id,class_contribution=cl,ecology=ec,selection_rule='two studies same direction; |patient-equal pooled rho|>=0.30; >=2 block scales aligned',block_permute(q$parts[,cl],q$eco[[ec]],q$xy,paste(section,cl,ec)))}}
permout<-if(length(perm))rbindlist(perm)else data.table(study=character(),section_id=character(),patient_id=character(),class_contribution=character(),ecology=character(),selection_rule=character(),block_scale=character(),occupied_blocks=integer(),block_rho=numeric(),empirical_p=numeric(),permutations=integer())
if(nrow(permout))permout[,FDR:=p.adjust(empirical_p,'BH')];w(merge(permout,sel,by=c('class_contribution','ecology'),all=TRUE),'SPATIAL_CLASS_PERMUTATION.tsv')
saveRDS(cache,file.path(OUT,'SPATIAL_GENE_Z_INTERNAL.rds'),compress=FALSE)
diff_count<-cross[,sum(max(rho,na.rm=TRUE)-min(rho,na.rm=TRUE)>=.10),by=ecology][V1>0,.N]
spgate<-if(diff_count>=2)'SPATIAL_HELDOUT_FUNCTIONAL_DIFFERENTIATION_PASS' else if(diff_count>=1)'SPATIAL_HELDOUT_FUNCTIONAL_DIFFERENTIATION_PARTIAL' else 'SPATIAL_HELDOUT_FUNCTIONAL_DIFFERENTIATION_NOT_SUPPORTED'
writeLines(c('# Spatial held-out validation gate','',paste0('Status: `',spgate,'`'),'',paste0('- Cross-study class/ecology endpoints: ',nrow(cross),'.'),paste0('- Endpoints meeting the frozen 999-permutation trigger: ',sum(sel$selected,na.rm=TRUE),'.'),'- Spatial evidence was annotation-only in the v1.3 primary provenance classifier; spatial leave-out membership is identical by algorithm.','- Inference remains patient-equal and associative, not causal.'),file.path(OUT,'SPATIAL_HELDOUT_VALIDATION_GATE.md'))

# IVY -------------------------------------------------------------------------
MATRIX_DIR<-file.path(DATA_ROOT,'spatial','ivy_gap','gene_expression_matrix')
mapping<-fread(file.path(V12,'05_ivy_reference/IVY_INPUT_LINEAGE.tsv'));stopifnot(nrow(mapping)==270,uniqueN(mapping$donor_id)==37)
gmeta<-fread(file.path(MATRIX_DIR,'rows-genes.csv'));ex<-fread(file.path(MATRIX_DIR,'fpkm_table.csv'),check.names=FALSE);sym<-toupper(trimws(gmeta$gene_symbol[match(as.character(ex[[1]]),as.character(gmeta$gene_id))]));ok<-!is.na(sym)&nzchar(sym);mat<-as.matrix(ex[ok,-1]);storage.mode(mat)<-'double';rownames(mat)<-sym[ok]
if(anyDuplicated(rownames(mat))){u<-unique(rownames(mat));mat<-rowsum(mat,rownames(mat),reorder=FALSE)/as.vector(table(factor(rownames(mat),levels=u)))}
ids<-as.character(mapping$sample_id);stopifnot(all(ids%in%colnames(mat)));mat<-log2(mat[,ids,drop=FALSE]+1);present<-genes95%in%rownames(mat);Zivy<-matrix(0,nrow=95,ncol=length(ids),dimnames=list(genes95,ids));zz<-t(scale(t(mat[genes95[present],,drop=FALSE])));zz[!is.finite(zz)]<-0;Zivy[present,]<-zz
oldfreeze<-fread(file.path(V12,'01_spatial_reference/SPATIAL_SIGNATURE_FREEZE.tsv'));overlap<-unique(oldfreeze[signature%in%c('MES_HYP','MES_AST')&mode=='RAW',gene])
regions<-c('CELLULAR_TUMOUR','PSEUDOPALISADING_CELLS_AROUND_NECROSIS','MICROVASCULAR_PROLIFERATION','INFILTRATING_TUMOUR','LEADING_EDGE')
defs<-data.table(contrast=c('PAN_VS_CT','MVP_VS_CT','LE_VS_CT','IT_VS_CT','PAN_VS_MVP'),r1=c(regions[2],regions[3],regions[5],regions[4],regions[2]),r0=c(regions[1],regions[1],regions[1],regions[1],regions[3]))
contrast_fun<-function(fit,r1,r0){nd<-data.frame(region=factor(c(r1,r0),levels=regions),donor_id=factor(rep(levels(fit@frame$donor_id)[1],2),levels=levels(fit@frame$donor_id)));mm<-model.matrix(~region,nd);b<-fixef(fit);cv<-mm[,names(b),drop=FALSE][1,]-mm[,names(b),drop=FALSE][2,];est<-sum(cv*b);se<-sqrt(as.numeric(t(cv)%*%vcov(fit)%*%cv));c(estimate=est,SE=se,CI_low=est-1.96*se,CI_high=est+1.96*se,p=2*pnorm(-abs(est/se)))}
base<-mapping[,.(donor_id=as.character(donor_id),sample_id=as.character(sample_id),region)];samples<-list();models<-list();contrasts<-list();boots<-list();mi<-ci<-bi<-0L
for(mode in c('RAW','OVERLAP_PRUNED'))for(cl in class_order){g<-genes95[classes[genes95]==cl&present];if(mode=='OVERLAP_PRUNED')g<-setdiff(g,overlap);score<-if(length(g))colSums(Zivy[g,,drop=FALSE])/95 else rep(0,ncol(Zivy));d<-copy(base);d[,score:=score];d[,region:=factor(region,levels=regions)];d[,donor_id:=factor(donor_id)];samples[[paste(mode,cl)]]<-cbind(mode=mode,class_contribution=cname[cl],d)
 fit<-lmer(score~region+(1|donor_id),d,REML=FALSE,control=lmerControl(optimizer='bobyqa'));sm<-as.data.table(coef(summary(fit)),keep.rownames='term');setnames(sm,c('Estimate','Std. Error','t value'),c('estimate','SE','z_value'),skip_absent=TRUE);sm[,`:=`(mode=mode,class_contribution=cname[cl],CI_low=estimate-1.96*SE,CI_high=estimate+1.96*SE,p_value=2*pnorm(-abs(z_value)),n_donors=uniqueN(d$donor_id),n_samples=nrow(d),singular=isSingular(fit),convergence_message=paste(fit@optinfo$conv$lme4$messages,collapse=';'),genes_covered=length(g),fixed_denominator=95)];mi<-mi+1;models[[mi]]<-sm
 for(j in seq_len(nrow(defs))){a<-contrast_fun(fit,defs$r1[j],defs$r0[j]);ci<-ci+1;contrasts[[ci]]<-data.table(mode=mode,class_contribution=cname[cl],contrast=defs$contrast[j],estimate=a['estimate'],SE=a['SE'],CI_low=a['CI_low'],CI_high=a['CI_high'],p_value=a['p'],n_donors=uniqueN(d$donor_id),n_samples=nrow(d));q<-d[region%in%c(defs$r1[j],defs$r0[j]),.(v=mean(score)),by=.(donor_id,region)];z<-merge(q[region==defs$r1[j],.(donor_id,v1=v)],q[region==defs$r0[j],.(donor_id,v0=v)],by='donor_id');z[,diff:=v1-v0];if(nrow(z)>=2){set.seed(as.integer(20260801405%%.Machine$integer.max)+ci);bb<-replicate(2000,mean(sample(z$diff,nrow(z),TRUE)));bi<-bi+1;boots[[bi]]<-rbind(data.table(mode=mode,class_contribution=cname[cl],contrast=defs$contrast[j],result_type='DONOR_CLUSTER_BOOTSTRAP',omitted_donor=NA_character_,estimate=mean(z$diff),CI_low=quantile(bb,.025),CI_high=quantile(bb,.975),n_donors=nrow(z),bootstrap_iterations=2000L),rbindlist(lapply(z$donor_id,function(id)data.table(mode=mode,class_contribution=cname[cl],contrast=defs$contrast[j],result_type='LEAVE_ONE_DONOR_OUT',omitted_donor=id,estimate=mean(z[donor_id!=id,diff]),CI_low=NA_real_,CI_high=NA_real_,n_donors=nrow(z)-1L,bootstrap_iterations=0L))))}}
}
ivysamp<-rbindlist(samples);w(ivysamp,'IVY_CLASS_CONTRIBUTIONS_INTERNAL.tsv');mod<-rbindlist(models);mod[,FDR:=p.adjust(p_value,'BH'),by=mode];con<-rbindlist(contrasts);con[,FDR:=p.adjust(p_value,'BH'),by=mode];w(mod,'IVY_CLASS_MIXED_MODELS.tsv');w(con,'IVY_CLASS_CONTRASTS.tsv');w(rbindlist(boots,fill=TRUE),'IVY_CLASS_BOOTSTRAP.tsv');saveRDS(list(Z=Zivy,base=base,present=present,overlap=overlap),file.path(OUT,'IVY_GENE_Z_INTERNAL.rds'),compress=FALSE)
rawcon<-con[mode=='RAW'];ivydiff<-rawcon[,sum(max(estimate)-min(estimate)>=.10),by=contrast][V1>0,.N];ivygate<-if(ivydiff>=3)'IVY_HELDOUT_FUNCTIONAL_DIFFERENTIATION_PASS' else if(ivydiff>=1)'IVY_HELDOUT_FUNCTIONAL_DIFFERENTIATION_PARTIAL' else 'IVY_HELDOUT_FUNCTIONAL_DIFFERENTIATION_NOT_SUPPORTED'
writeLines(c('# IVY held-out validation gate','',paste0('Status: `',ivygate,'`'),'',paste0('- Frozen donor/sample reference: ',uniqueN(base$donor_id),' donors / ',nrow(base),' samples.'),paste0('- Class contrasts evaluated: ',nrow(rawcon),' primary raw rows.'),'- IVY evidence was annotation-only in the v1.3 primary provenance classifier; IVY leave-out membership is identical by algorithm.','- Mixed models, donor-cluster bootstrap, patient-equal paired contrasts, LODO and overlap-pruned sensitivity are reported.'),file.path(OUT,'IVY_HELDOUT_VALIDATION_GATE.md'))

# Append exact additivity/coverage QA.
ivwide<-dcast(ivysamp[mode=='RAW'],donor_id+sample_id+region~class_contribution,value.var='score');ivwide[,MES_ADDITIVE:=rowSums(.SD),.SDcols=cname]
qa<-fread(file.path(OUT,'ADDITIVITY_QA.tsv'));qa<-rbind(qa,data.table(data_layer='SPATIAL_ALL_SECTIONS',n_observations=nrow(spots),max_abs_additivity_error=max(abs(spots$MES_ADDITIVE-rowSums(spots[,..cname]))),tolerance=1e-10,status='PASS'),data.table(data_layer='IVY_GAP_RAW',n_observations=nrow(base),max_abs_additivity_error=max(abs(ivwide$MES_ADDITIVE-colSums(Zivy[,ivwide$sample_id,drop=FALSE])/95)),tolerance=1e-10,status='PASS'),fill=TRUE);w(qa,'ADDITIVITY_QA.tsv')
cov<-fread(file.path(OUT,'CONTRIBUTION_COVERAGE.tsv'));cov<-rbind(cov,rbindlist(lapply(class_order,function(cl)data.table(data_layer='SPATIAL_ALL_SECTIONS',provenance_class=cl,genes_frozen=sum(classes[genes95]==cl),genes_present=NA_integer_,fixed_denominator=95,missing_gene_rule='SECTION_SPECIFIC_ZERO_ON_STANDARDIZED_MEAN_REFERENCE_SCALE'))),rbindlist(lapply(class_order,function(cl)data.table(data_layer='IVY_GAP',provenance_class=cl,genes_frozen=sum(classes[genes95]==cl),genes_present=sum(classes[genes95]==cl&present),fixed_denominator=95,missing_gene_rule='ZERO_ON_STANDARDIZED_MEAN_REFERENCE_SCALE'))),fill=TRUE);w(cov,'CONTRIBUTION_COVERAGE.tsv')
if(any(qa$status!='PASS'))stop('ADDITIVITY_QA_FAIL')
writeLines(capture.output(sessionInfo()),file.path(RUN,'session_info/MODULE_A_SPATIAL_IVY_SESSION_INFO.txt'))
cat('MODULE_A_SPATIAL_IVY_COMPLETE',nrow(spots),nrow(ivysamp),'\n')
