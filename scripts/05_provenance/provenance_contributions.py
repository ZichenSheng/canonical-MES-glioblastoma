#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import platform
import re
import os
from pathlib import Path

DATA_ROOT = Path(os.environ.get("GBM_MES_DATA_ROOT", "data"))
RESULT_ROOT = Path(os.environ.get("GBM_MES_OUTPUT_ROOT", "results"))
REPO_ROOT = Path(os.environ.get("GBM_MES_REPO_ROOT", Path(__file__).resolve().parents[2]))

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.optimize import minimize_scalar

RUN = RESULT_ROOT / 'provenance'
V13 = RESULT_ROOT / 'interpretability'
V12 = RESULT_ROOT / 'spatial'
V11 = RESULT_ROOT / 'scoring'
V1 = DATA_ROOT / 'prepared' / 'core'
DATA = DATA_ROOT / 'prepared' / 'single_cell_reference'
SIG = REPO_ROOT / 'resources' / 'signatures' / 'signature_gene_sets.tsv'
OUT = RUN / 'contributions'
CLASSES = ['MALIGNANT_DOMINANT', 'ECOLOGICAL_DOMINANT', 'SHARED_MIXED_ORIGIN', 'LOW_INFORMATION_OR_UNSTABLE']
CNAMES = {'MALIGNANT_DOMINANT':'C_MALIGNANT', 'ECOLOGICAL_DOMINANT':'C_ECOLOGICAL',
          'SHARED_MIXED_ORIGIN':'C_SHARED', 'LOW_INFORMATION_OR_UNSTABLE':'C_UNSTABLE'}
SEEDS = json.loads((RUN/'00_governance/V1_4_RANDOM_SEEDS.json').read_text())


def wt(df: pd.DataFrame, name: str) -> None:
    df.to_csv(OUT/name, sep='\t', index=False, na_rep='NA')


def zcols(expr: pd.DataFrame, genes: list[str]) -> np.ndarray:
    out = np.zeros((expr.shape[1], len(genes)), dtype=np.float64)
    for j, gene in enumerate(genes):
        if gene not in expr.index:
            continue
        x = expr.loc[gene].to_numpy(float)
        sd = np.nanstd(x, ddof=1)
        if np.isfinite(sd) and sd > 0:
            out[:, j] = (x - np.nanmean(x)) / sd
    return out


def score_parts(z: np.ndarray, genes: list[str], cmap: dict[str,str], present: set[str], pruned: set[str] | None = None) -> pd.DataFrame:
    pruned = pruned or set()
    out = {}
    for cl in CLASSES:
        idx = [j for j,g in enumerate(genes) if cmap[g] == cl and g in present and g not in pruned]
        out[CNAMES[cl]] = z[:,idx].sum(axis=1)/95.0 if idx else np.zeros(z.shape[0])
    d = pd.DataFrame(out)
    d['MES_ADDITIVE'] = d[list(CNAMES.values())].sum(axis=1)
    return d


def boot_spearman(x: np.ndarray, y: np.ndarray, seed: int, B: int=2000) -> tuple:
    ok = np.isfinite(x)&np.isfinite(y); x=x[ok]; y=y[ok]; n=len(x)
    if n < 4: return (n,np.nan,np.nan,np.nan,np.nan)
    rho = stats.spearmanr(x,y).statistic
    rng=np.random.default_rng(seed); vals=np.empty(B)
    for b in range(B):
        ii=rng.integers(0,n,n); vals[b]=stats.spearmanr(x[ii],y[ii]).statistic
    return n,rho,*np.nanquantile(vals,[.025,.975]),np.nanstd(vals,ddof=1)


def reml_meta(yi: np.ndarray, vi: np.ndarray) -> dict:
    ok=np.isfinite(yi)&np.isfinite(vi)&(vi>0); y=yi[ok]; v=vi[ok]; k=len(y)
    if k == 0: return dict(k=0,estimate=np.nan,ci_low=np.nan,ci_high=np.nan,prediction_low=np.nan,prediction_high=np.nan,tau2=np.nan,I2=np.nan)
    if k == 1:
        se=math.sqrt(v[0]); return dict(k=1,estimate=y[0],ci_low=y[0]-1.96*se,ci_high=y[0]+1.96*se,prediction_low=np.nan,prediction_high=np.nan,tau2=0,I2=0)
    def obj(t):
        w=1/(v+t); mu=np.sum(w*y)/np.sum(w)
        return .5*(np.sum(np.log(v+t))+np.log(np.sum(w))+np.sum(w*(y-mu)**2))
    opt=minimize_scalar(obj,bounds=(0,max(10,float(np.var(y)*10+1))),method='bounded')
    tau=max(0,float(opt.x)); w=1/(v+tau); mu=float(np.sum(w*y)/np.sum(w)); se=math.sqrt(1/np.sum(w))
    q=float(np.sum((1/v)*(y-np.sum(y/v)/np.sum(1/v))**2)); i2=max(0,(q-(k-1))/q*100) if q>0 else 0
    return dict(k=k,estimate=mu,ci_low=mu-1.96*se,ci_high=mu+1.96*se,
                prediction_low=mu-1.96*math.sqrt(tau+se*se),prediction_high=mu+1.96*math.sqrt(tau+se*se),tau2=tau,I2=i2)


def load_sets() -> tuple[list[str],dict[str,str],set[str]]:
    p=pd.read_csv(RUN/'01_provenance_freeze/PROVENANCE_CLASS_FREEZE.tsv',sep='\t')
    gcol='gene' if 'gene' in p.columns else 'harmonized_symbol'; ccol='primary_origin' if 'primary_origin' in p.columns else 'origin_class'
    p[gcol]=p[gcol].astype(str).str.upper(); genes=p[gcol].tolist(); cmap=dict(zip(p[gcol],p[ccol]))
    reg=pd.read_csv(SIG,sep='\t'); reg['gene']=reg.gene_symbol_clean.astype(str).str.strip().str.upper()
    ctx=set(reg.loc[reg.signature_id.isin(['HALLMARK_HYPOXIA','HALLMARK_EMT','NABA_CORE_MATRISOME','MYELOID_CORE_IDENTITY']),'gene'])
    return genes,cmap,ctx


def read_bulk_matrix(path: Path, ids: list[str], log2_transform: bool) -> pd.DataFrame:
    hdr=pd.read_csv(path,sep='\t',nrows=0,compression='infer').columns.tolist()
    if not set(ids).issubset(hdr): raise RuntimeError(f'missing frozen bulk IDs in {path}')
    d=pd.read_csv(path,sep='\t',usecols=[hdr[0]]+ids,compression='infer')
    d.iloc[:,0]=d.iloc[:,0].astype(str).str.strip().str.upper()
    x=d.groupby(d.columns[0],sort=False)[ids].mean()
    if log2_transform: x=np.log2(x.astype(float)+1)
    return x.astype(float)


def strict_bulk(genes: list[str], cmap: dict[str,str], ctx_overlap: set[str]) -> tuple[pd.DataFrame,dict[str,np.ndarray],dict[str,list[str]]]:
    mem=DATA_ROOT/'prepared/cohorts'
    specs={
      'CGGA325':(DATA_ROOT/'bulk/cgga325_expression.tsv','membership_cgga325.tsv',True),
      'CGGA693':(DATA_ROOT/'bulk/cgga693_expression.tsv','membership_cgga693.tsv',True),
      'TCGA':(DATA_ROOT/'bulk/tcga_gbm_expression.tsv.gz','membership_tcga.tsv',False)}
    prior=pd.read_csv(V1/'05_bulk/BULK_PATIENT_SCORES.tsv',sep='\t')
    tcga=pd.read_csv(V11/'06_bulk_integration/TCGA_BULK_PATIENT_SCORES.tsv',sep='\t');tcga['cohort']='TCGA';prior=pd.concat([prior,tcga],ignore_index=True,sort=False)
    bp=pd.read_csv(V13/'01_state_controls/BAYESPRISM_STATE_PATIENT_SCORES_INTERNAL.tsv',sep='\t')
    rows=[]; assoc=[]; matrices={}; ids_by={}; qa=[]; coverage=[]
    for ci,(co,(path,mfile,do_log)) in enumerate(specs.items()):
        md=pd.read_csv(mem/mfile,sep='\t'); idcol=next(c for c in ['patient_id','sample_id','entity_id'] if c in md.columns); ids=md[idcol].astype(str).tolist()
        if co=='TCGA':
            tm=pd.read_csv(V11/'05_tcga_strict/TCGA_STRICT_COHORT.tsv',sep='\t')
            tm=tm[tm.strict_inclusion.astype(bool) & tm.patient_id.astype(str).isin(ids)]
            if len(tm)!=len(ids): raise RuntimeError('TCGA frozen patient-to-sample mapping mismatch')
            ids=tm.sample_id.astype(str).tolist()
        expr=read_bulk_matrix(path,ids,do_log); present=set(genes)&set(expr.index); z=zcols(expr,genes); matrices[co]=z; ids_by[co]=ids
        parts=score_parts(z,genes,cmap,present); pr=score_parts(z,genes,cmap,present,ctx_overlap)
        cov=prior[(prior.cohort==co)&prior.patient_id.astype(str).isin(ids)].copy(); cov['patient_id']=cov.patient_id.astype(str); cov=cov.set_index('patient_id').reindex(ids).reset_index()
        if co!='TCGA':
            b=bp[(bp.cohort==co)&bp.patient_id.astype(str).isin(ids)][['patient_id','myeloid_fraction','vascular_stromal_fraction','posterior_median_cv']].copy();b['patient_id']=b.patient_id.astype(str);cov=cov.merge(b,on='patient_id',how='left')
        else:
            cov['myeloid_fraction']=np.nan;cov['vascular_stromal_fraction']=np.nan;cov['posterior_median_cv']=np.nan
        d=pd.concat([pd.DataFrame({'cohort':co,'patient_id':ids}),parts],axis=1)
        for c in CNAMES.values(): d[c+'_PRUNED']=pr[c]
        for c in ['context','hypoxia','EMT','matrix','myeloid','myeloid_fraction','vascular_stromal_fraction','posterior_median_cv']:
            d[c]=cov[c].to_numpy() if c in cov else np.nan
        top=d.MES_ADDITIVE>=d.MES_ADDITIVE.quantile(.75); denom=d.loc[top,list(CNAMES.values())].abs().mean().sum()
        for c in CNAMES.values():
            n,r,lo,hi,se=boot_spearman(d[c].to_numpy(),d.MES_ADDITIVE.to_numpy(),SEEDS['strict_bulk_bootstrap']+ci*100+list(CNAMES.values()).index(c))
            assoc.append(dict(cohort=co,analysis='RAW',class_contribution=c,endpoint='TOTAL_MES',n=n,rho=r,ci_low=lo,ci_high=hi,bootstrap_se=se,bootstrap_replicates=2000))
            n,r,lo,hi,se=boot_spearman(d[c].to_numpy(),d.context.to_numpy(),SEEDS['strict_bulk_bootstrap']+500+ci*100+list(CNAMES.values()).index(c))
            assoc.append(dict(cohort=co,analysis='RAW',class_contribution=c,endpoint='CONTEXT',n=n,rho=r,ci_low=lo,ci_high=hi,bootstrap_se=se,bootstrap_replicates=2000))
            assoc.append(dict(cohort=co,analysis='RAW',class_contribution=c,endpoint='MES_TOP_QUARTILE_SIGNED_MEAN',n=int(top.sum()),rho=float(d.loc[top,c].mean()),ci_low=np.nan,ci_high=np.nan,bootstrap_se=np.nan,bootstrap_replicates=0))
            assoc.append(dict(cohort=co,analysis='RAW',class_contribution=c,endpoint='MES_TOP_QUARTILE_ABSOLUTE_RELATIVE_COMPOSITION',n=int(top.sum()),rho=float(d.loc[top,c].abs().mean()/denom) if denom else np.nan,ci_low=np.nan,ci_high=np.nan,bootstrap_se=np.nan,bootstrap_replicates=0))
            n,r,lo,hi,se=boot_spearman(d[c+'_PRUNED'].to_numpy(),d.context.to_numpy(),SEEDS['strict_bulk_bootstrap']+1000+ci*100+list(CNAMES.values()).index(c))
            assoc.append(dict(cohort=co,analysis='OVERLAP_PRUNED_FIXED_DENOMINATOR',class_contribution=c,endpoint='CONTEXT',n=n,rho=r,ci_low=lo,ci_high=hi,bootstrap_se=se,bootstrap_replicates=2000))
        err=np.max(np.abs(d.MES_ADDITIVE-d[list(CNAMES.values())].sum(axis=1)))
        qa.append(dict(data_layer=f'STRICT_BULK_{co}',n_observations=len(d),max_abs_additivity_error=err,tolerance=1e-10,status='PASS' if err<1e-10 else 'FAIL'))
        for cl in CLASSES:
            coverage.append(dict(data_layer=f'STRICT_BULK_{co}',provenance_class=cl,genes_frozen=sum(cmap[g]==cl for g in genes),genes_present=sum(cmap[g]==cl and g in present for g in genes),fixed_denominator=95,missing_gene_rule='ZERO_ON_STANDARDIZED_MEAN_REFERENCE_SCALE'))
        rows.append(d)
    full=pd.concat(rows,ignore_index=True); ass=pd.DataFrame(assoc)
    metas=[]
    for (analysis,c,ep),q in ass[ass.endpoint.isin(['TOTAL_MES','CONTEXT'])].groupby(['analysis','class_contribution','endpoint']):
        yi=np.arctanh(np.clip(q.rho.to_numpy(float),-.999999,.999999)); vi=1/(q.n.to_numpy(float)-3); mm=reml_meta(yi,vi)
        for key in ['estimate','ci_low','ci_high','prediction_low','prediction_high']: mm[key]=np.tanh(mm[key]) if np.isfinite(mm[key]) else np.nan
        mm.update(analysis=analysis,class_contribution=c,endpoint=ep,direction_consistent=len(set(np.sign(q.rho)))==1,cohort_effects=';'.join(f'{a}:{b:.6g}' for a,b in zip(q.cohort,q.rho)))
        metas.append(mm)
    wt(full,'STRICT_BULK_CLASS_CONTRIBUTIONS.tsv');wt(ass,'STRICT_BULK_CLASS_CONTEXT_ASSOCIATIONS.tsv');wt(pd.DataFrame(metas),'STRICT_BULK_CLASS_META.tsv')
    np.savez_compressed(OUT/'STRICT_BULK_GENE_Z_INTERNAL.npz',genes=np.array(genes),**{f'Z_{k}':v for k,v in matrices.items()},**{f'IDS_{k}':np.array(v) for k,v in ids_by.items()})
    return pd.DataFrame(qa),matrices,ids_by,pd.DataFrame(coverage)


def fit_effect(d: pd.DataFrame, exp: str) -> dict:
    y=d.score.to_numpy(float); mes=d.malignant_mes_fraction.to_numpy(float); comp=d.composition_fraction.to_numpy(float)
    pat=pd.get_dummies(d.patient_key,drop_first=True,dtype=float).to_numpy()
    if exp=='A_FIXED_MALIGNANT': terms=np.column_stack([comp]); names=['composition_slope']
    elif exp=='B_FIXED_COMPOSITION': terms=np.column_stack([mes]); names=['malignant_state_slope']
    else: terms=np.column_stack([mes,comp,mes*comp]); names=['malignant_state_slope','composition_slope','interaction']
    X=np.column_stack([np.ones(len(d)),terms,pat]);keep=np.isfinite(y)&np.isfinite(X).all(axis=1);X=X[keep];y=y[keep]
    beta=np.linalg.lstsq(X,y,rcond=None)[0]; resid=y-X@beta; df=max(1,len(y)-np.linalg.matrix_rank(X)); s2=np.sum(resid**2)/df; cov=s2*np.linalg.pinv(X.T@X)
    out={'n':len(y),'patients':d.patient_key.nunique(),'studies':d.study.nunique(),'model_R2':1-np.sum(resid**2)/np.sum((y-y.mean())**2)}
    for j,nm in enumerate(names,1): out[nm]=beta[j];out[nm+'_ci_low']=beta[j]-1.96*math.sqrt(max(0,cov[j,j]));out[nm+'_ci_high']=beta[j]+1.96*math.sqrt(max(0,cov[j,j]))
    psl=[];ssl=[]
    for _,q in d.groupby('patient_key'):
        xx=np.column_stack([np.ones(len(q)), q.composition_fraction]) if exp=='A_FIXED_MALIGNANT' else np.column_stack([np.ones(len(q)),q.malignant_mes_fraction]) if exp=='B_FIXED_COMPOSITION' else np.column_stack([np.ones(len(q)),q.malignant_mes_fraction,q.composition_fraction,q.malignant_mes_fraction*q.composition_fraction])
        if len(q)>xx.shape[1] and np.linalg.matrix_rank(xx)==xx.shape[1]: psl.append(np.linalg.lstsq(xx,q.score,rcond=None)[0][1:])
    for _,q in d.groupby('study'):
        xx=np.column_stack([np.ones(len(q)), q.composition_fraction]) if exp=='A_FIXED_MALIGNANT' else np.column_stack([np.ones(len(q)),q.malignant_mes_fraction]) if exp=='B_FIXED_COMPOSITION' else np.column_stack([np.ones(len(q)),q.malignant_mes_fraction,q.composition_fraction,q.malignant_mes_fraction*q.composition_fraction])
        if len(q)>xx.shape[1] and np.linalg.matrix_rank(xx)==xx.shape[1]: ssl.append(np.linalg.lstsq(xx,q.score,rcond=None)[0][1:])
    out['patient_slope_heterogeneity_sd']=float(np.nanmedian(np.nanstd(np.asarray(psl),axis=0,ddof=1))) if len(psl)>1 else np.nan
    out['study_slope_heterogeneity_sd']=float(np.nanmedian(np.nanstd(np.asarray(ssl),axis=0,ddof=1))) if len(ssl)>1 else np.nan
    return out


def pseudobulk(genes: list[str], cmap: dict[str,str]) -> tuple[pd.DataFrame,pd.DataFrame]:
    meta=pd.read_csv(DATA/'FINAL_GBM_BAYESPRISM_REFERENCE_METADATA.tsv',sep='\t')
    reg=pd.read_csv(SIG,sep='\t');reg['gene']=reg.gene_symbol_clean.astype(str).str.strip().str.upper()
    formal_ids=['NEFTEL_MES_LIKE','NEFTEL_MES1','NEFTEL_MES2','NEFTEL_AC_LIKE','NEFTEL_OPC_LIKE','NEFTEL_NPC_LIKE','PROLIFERATION_CELL_CYCLE']
    wanted=set(reg.loc[reg.signature_id.isin(formal_ids),'gene'].dropna())
    wanted.update(reg.loc[reg.signature_id.isin(['NEFTEL_MES_HYP','NEFTEL_MES_AST']),'gene'].dropna())
    with h5py.File(DATA/'REF_B_MULTISTUDY_SCRNA_COUNTS.h5','r') as h:
        ids=h['cell_id'].asstr()[:]; allg=np.char.upper(h['gene_symbol'].asstr()[:].astype(str));grp=h['counts_csr_cells_by_genes'];full=sparse.csr_matrix((grp['data'][:],grp['indices'][:],grp['indptr'][:]),shape=(len(ids),len(allg)))
    if not np.array_equal(ids,meta.cell_id.astype(str).to_numpy()): raise RuntimeError('pseudo-bulk reference identity mismatch')
    fmap={g:i for i,g in enumerate(allg)}; wanted=sorted(wanted & set(allg));present=[g for g in genes if g in fmap]; x=full[:,[fmap[g] for g in wanted]].tocsr();del full
    meta['patient_key']=meta.author.astype(str)+'::'+meta.donor_id.astype(str);meta['row']=np.arange(len(meta))
    excluded={'Neftel2019','Bhaduri2020','Couturier2020','Johnson2020','Richards2021','Yu2020','Yuan2018'};pool=meta[~meta.author.isin(excluded)].copy();pools={}
    for patient,d in pool.groupby('patient_key'):
        pools[patient]={'study':d.author.iloc[0],'malignant':d.loc[d.reference_cell_type.eq('malignant'),'row'].to_numpy(int),'MES':d.loc[d.reference_cell_type.eq('malignant')&d.malignant_state.eq('MES-like'),'row'].to_numpy(int),'nonMES':d.loc[d.reference_cell_type.eq('malignant')&~d.malignant_state.eq('MES-like'),'row'].to_numpy(int),'myeloid':d.loc[d.reference_cell_type.eq('myeloid_macrophage'),'row'].to_numpy(int),'endothelial':d.loc[d.reference_cell_type.eq('endothelial'),'row'].to_numpy(int),'pericyte':d.loc[d.reference_cell_type.eq('pericyte'),'row'].to_numpy(int),'lymphoid':d.loc[d.reference_cell_type.isin(['t_cell','nk_cell','b_plasma']),'row'].to_numpy(int),'oligodendrocyte_other':d.loc[d.reference_cell_type.isin(['oligodendrocyte','opc_nonmalignant','mast_cell']),'row'].to_numpy(int),'all_nonmalignant':d.loc[~d.reference_cell_type.eq('malignant'),'row'].to_numpy(int)}
    rng=np.random.default_rng(SEEDS['pseudobulk_regeneration']);grid=np.array([0,.05,.10,.20,.30,.40]);reps=10
    def samp(a,n): return np.empty(0,dtype=int) if n<=0 else rng.choice(a,n,replace=n>len(a))
    designs=[];vectors=[];sid=0
    for comp in ['myeloid','endothelial','pericyte','lymphoid','oligodendrocyte_other']:
      for patient,p in pools.items():
       if len(p['malignant'])<20 or len(p[comp])<10: continue
       fixed=samp(p['malignant'],100);mes0=np.mean(meta.iloc[fixed].malignant_state.eq('MES-like'))
       for frac in grid:
        ncomp=int(round(100*frac/(1-frac)))
        for rep in range(reps):
         rr=np.r_[fixed,samp(p[comp],ncomp)];sid+=1;vectors.append(np.asarray(x[rr].sum(axis=0)).ravel());designs.append(dict(pseudobulk_id=f'PB{sid:06d}',experiment='A_FIXED_MALIGNANT',compartment=comp,study=p['study'],patient_key=patient,replicate=rep+1,target_fraction=frac,malignant_mes_fraction=mes0,composition_fraction=ncomp/len(rr),n_malignant=100,n_nonmalignant=ncomp,integer_count_source='REF_B raw-count H5; source studies exclude discovery and current validation'))
    for patient,p in pools.items():
      if len(p['MES'])<10 or len(p['nonMES'])<10 or len(p['all_nonmalignant'])<10:continue
      fixed=samp(p['all_nonmalignant'],25)
      for frac in grid:
       nmes=int(round(100*frac));nnon=100-nmes
       for rep in range(reps):
        rr=np.r_[samp(p['MES'],nmes),samp(p['nonMES'],nnon),fixed];sid+=1;vectors.append(np.asarray(x[rr].sum(axis=0)).ravel());designs.append(dict(pseudobulk_id=f'PB{sid:06d}',experiment='B_FIXED_COMPOSITION',compartment='all_nonmalignant',study=p['study'],patient_key=patient,replicate=rep+1,target_fraction=frac,malignant_mes_fraction=nmes/100,composition_fraction=.2,n_malignant=100,n_nonmalignant=25,integer_count_source='REF_B raw-count H5; source studies exclude discovery and current validation'))
    for comp in ['myeloid','endothelial','pericyte']:
     for patient,p in pools.items():
      if len(p['MES'])<10 or len(p['nonMES'])<10 or len(p[comp])<10:continue
      for mf in grid:
       nmes=int(round(100*mf));nnon=100-nmes
       for cf in grid:
        ncomp=int(round(100*cf/(1-cf)))
        for rep in range(reps):
         rr=np.r_[samp(p['MES'],nmes),samp(p['nonMES'],nnon),samp(p[comp],ncomp)];sid+=1;vectors.append(np.asarray(x[rr].sum(axis=0)).ravel());designs.append(dict(pseudobulk_id=f'PB{sid:06d}',experiment='C_TWO_DIMENSIONAL',compartment=comp,study=p['study'],patient_key=patient,replicate=rep+1,target_fraction=f'MES={mf};COMP={cf}',malignant_mes_fraction=mf,composition_fraction=ncomp/len(rr),n_malignant=100,n_nonmalignant=ncomp,integer_count_source='REF_B raw-count H5; source studies exclude discovery and current validation'))
    design=pd.DataFrame(designs);formal=pd.read_csv(V11/'04_pseudobulk_benchmark/PSEUDOBULK_DESIGN.tsv',sep='\t')
    design_ok=len(design)==len(formal) and design.pseudobulk_id.equals(formal.pseudobulk_id) and design[['experiment','compartment','study','patient_key']].equals(formal[['experiment','compartment','study','patient_key']]) and np.allclose(design[['malignant_mes_fraction','composition_fraction']].astype(float),formal[['malignant_mes_fraction','composition_fraction']].astype(float),atol=1e-12)
    if not design_ok: raise RuntimeError('exact formal pseudo-bulk design regeneration failed')
    counts=np.vstack(vectors);lib=counts.sum(axis=1);logcpm=np.log2((counts+.5)/(lib[:,None]+1)*1e6);Z=np.zeros((len(design),95));present_idx=[genes.index(g) for g in present];wanted_idx=[wanted.index(g) for g in present]
    for (_, _),ii0 in design.groupby(['experiment','compartment']).groups.items():
        ii=np.asarray(list(ii0));a=logcpm[ii];mu=a.mean(axis=0);sd=a.std(axis=0,ddof=1);sd[sd==0]=np.nan;zz=(a-mu)/sd;zz[~np.isfinite(zz)]=0;Z[np.ix_(ii,present_idx)]=zz[:,wanted_idx]
    original=Z[:,present_idx].mean(axis=1);resp=pd.read_csv(V11/'04_pseudobulk_benchmark/PSEUDOBULK_SCORE_RESPONSES.tsv',sep='\t');resp=resp[resp.score_name.eq('canonical_MES')].set_index('pseudobulk_id').reindex(design.pseudobulk_id)
    # The frozen response TSV is serialized to seven decimal places; 5e-7 is
    # the strict rounding envelope. Additivity itself remains tested at 1e-10.
    max_repro=float(np.max(np.abs(original-resp.score.to_numpy())));repro_ok=max_repro<5e-7
    parts=score_parts(Z,genes,cmap,set(present));base=pd.concat([design.reset_index(drop=True),parts],axis=1)
    results=[]
    for cls in CNAMES.values():
      for (exp,comp),q in base.groupby(['experiment','compartment']):
        dd=q.copy();dd['score']=dd[cls];r=fit_effect(dd,exp);r.update(experiment=exp,compartment=comp,class_contribution=cls,classification='FROZEN_V1_3',scale='FIXED_95_DENOMINATOR');results.append(r)
    frozen=pd.DataFrame(results);leave=frozen.copy();leave['classification']='PSEUDOBULK_LEAVEOUT_V1_3_ALGORITHM';leave['leaveout_effect_on_primary_class']='NONE_PSEUDOBULK_WAS_ANNOTATION_ONLY'
    wt(frozen[frozen.experiment!='C_TWO_DIMENSIONAL'],'PSEUDOBULK_FROZEN_CLASS_RESULTS.tsv');wt(leave[leave.experiment!='C_TWO_DIMENSIONAL'],'PSEUDOBULK_LEAVEOUT_CLASS_RESULTS.tsv');wt(pd.concat([frozen[frozen.experiment.eq('C_TWO_DIMENSIONAL')],leave[leave.experiment.eq('C_TWO_DIMENSIONAL')]],ignore_index=True),'PSEUDOBULK_CLASS_INTERACTIONS.tsv')
    np.savez_compressed(OUT/'PSEUDOBULK_GENE_Z_INTERNAL.npz',genes=np.array(genes),Z=Z,pseudobulk_id=design.pseudobulk_id.to_numpy())
    state=frozen[frozen.experiment.eq('B_FIXED_COMPOSITION')].set_index('class_contribution').malignant_state_slope.abs();comp=frozen[frozen.experiment.eq('A_FIXED_MALIGNANT')].groupby('class_contribution').composition_slope.apply(lambda x:np.median(np.abs(x)))
    state_ok=np.nanmedian([state.get('C_MALIGNANT',np.nan),state.get('C_SHARED',np.nan)])>state.get('C_ECOLOGICAL',np.inf);comp_ok=np.nanmedian([comp.get('C_ECOLOGICAL',np.nan),comp.get('C_SHARED',np.nan)])>comp.get('C_MALIGNANT',np.inf)
    gate='PSEUDOBULK_HELDOUT_VALIDATION_PASS' if state_ok and comp_ok else 'PSEUDOBULK_HELDOUT_VALIDATION_PARTIAL' if state_ok or comp_ok else 'PSEUDOBULK_HELDOUT_VALIDATION_NOT_SUPPORTED'
    (OUT/'PSEUDOBULK_HELDOUT_VALIDATION_GATE.md').write_text(f'# Pseudo-bulk held-out validation gate\n\nStatus: `{gate}`\n\n- Exact frozen design regeneration: {design_ok}.\n- Canonical detected-gene response max absolute reproduction error: {max_repro:.3g} (PASS={repro_ok}).\n- Malignant/shared state-response ordering criterion: {state_ok}.\n- Ecological/shared composition-response ordering criterion: {comp_ok}.\n- Pseudo-bulk leave-out classes equal frozen primary classes because pseudo-bulk was annotation-only in the v1.3 classifier.\n')
    qa=pd.DataFrame([dict(data_layer='PSEUDOBULK_FORMAL',n_observations=len(base),max_abs_additivity_error=float(np.max(np.abs(base.MES_ADDITIVE-base[list(CNAMES.values())].sum(axis=1)))),tolerance=1e-10,status='PASS'),dict(data_layer='PSEUDOBULK_CANONICAL_RESPONSE_REPRODUCTION',n_observations=len(base),max_abs_additivity_error=max_repro,tolerance=5e-7,status='PASS' if repro_ok else 'FAIL')])
    cov=pd.DataFrame([dict(data_layer='PSEUDOBULK_FORMAL',provenance_class=cl,genes_frozen=sum(cmap[g]==cl for g in genes),genes_present=sum(cmap[g]==cl and g in present for g in genes),fixed_denominator=95,missing_gene_rule='ZERO_ON_STANDARDIZED_MEAN_REFERENCE_SCALE') for cl in CLASSES])
    return qa,cov


def discordance(bulk: pd.DataFrame) -> None:
    dis=pd.read_csv(V13/'02_patient_discordance/PATIENT_LEVEL_DISCORDANCE.tsv',sep='\t');d=dis.merge(bulk[['cohort','patient_id','MES_ADDITIVE']+list(CNAMES.values())],on=['cohort','patient_id'],how='left',validate='one_to_one')
    d['B_EM']=d.C_ECOLOGICAL-d.C_MALIGNANT;d['B_SM']=d.C_SHARED-d.C_MALIGNANT
    wt(d[['cohort','patient_id','MES_ADDITIVE']+list(CNAMES.values())+['D_residual','D_difference','quadrant','B_EM','B_SM']],'DISCORDANCE_CLASS_CONTRIBUTIONS.tsv')
    endpoints=list(CNAMES.values())+['B_EM','B_SM'];rows=[];rng=np.random.default_rng(SEEDS['discordance_bootstrap'])
    for co,q in d.groupby('cohort'):
      for ep in endpoints:
        x=q[ep].to_numpy(float);y=q.D_residual.to_numpy(float);ok=np.isfinite(x)&np.isfinite(y);x=x[ok];y=y[ok];X=np.column_stack([np.ones(len(x)),x]);beta=np.linalg.lstsq(X,y,rcond=None)[0][1];rho=stats.spearmanr(x,y).statistic;bs=[]
        for _ in range(2000):
          ii=rng.integers(0,len(x),len(x));bs.append(np.linalg.lstsq(np.column_stack([np.ones(len(x)),x[ii]]),y[ii],rcond=None)[0][1])
        lo,hi=np.quantile(bs,[.025,.975]);rows.append(dict(cohort=co,predictor=ep,outcome='D_residual_PRIMARY',n=len(x),spearman_rho=rho,robust_slope=beta,ci_low=lo,ci_high=hi,p_value=2*min(np.mean(np.asarray(bs)<=0),np.mean(np.asarray(bs)>=0)),bootstrap_replicates=2000,model='UNIVARIABLE_PATIENT_BOOTSTRAP'))
    models=pd.DataFrame(rows);models['FDR']=stats.false_discovery_control(models.p_value.fillna(1).to_numpy(),method='bh');wt(models,'DISCORDANCE_BALANCE_MODELS.tsv')
    comps=[]
    for co,q in d.groupby('cohort'):
      a=q[q.quadrant.eq('bulk-high / malignant-low')];b=q[q.quadrant.eq('bulk-high / malignant-high')]
      for ep in list(CNAMES.values())+['B_EM','B_SM']:
        est=a[ep].mean()-b[ep].mean();p=stats.mannwhitneyu(a[ep],b[ep],alternative='two-sided').pvalue if len(a) and len(b) else np.nan;comps.append(dict(cohort=co,contrast='bulk-high_malignant-low_MINUS_bulk-high_malignant-high',endpoint=ep,n_group1=len(a),n_group0=len(b),mean_group1=a[ep].mean(),mean_group0=b[ep].mean(),mean_difference=est,p_value=p))
    comp=pd.DataFrame(comps);comp['FDR']=stats.false_discovery_control(comp.p_value.fillna(1).to_numpy(),method='bh');wt(comp,'DISCORDANCE_QUADRANT_COMPARISONS.tsv')
    metas=[]
    for ep,q in models.groupby('predictor'):
        se=(q.ci_high-q.ci_low)/(2*1.96);mm=reml_meta(q.robust_slope.to_numpy(),se.to_numpy()**2);mm.update(predictor=ep,outcome='D_residual_PRIMARY',direction_consistent=len(set(np.sign(q.robust_slope)))==1);metas.append(mm)
    wt(pd.DataFrame(metas),'DISCORDANCE_CLASS_META.tsv')


def main() -> None:
    genes,cmap,ctx=load_sets()
    qa1,_,_,cov1=strict_bulk(genes,cmap,ctx)
    bulk=pd.read_csv(OUT/'STRICT_BULK_CLASS_CONTRIBUTIONS.tsv',sep='\t')
    discordance(bulk[bulk.cohort.isin(['CGGA325','CGGA693'])])
    qa2,cov2=pseudobulk(genes,cmap)
    qa=pd.concat([qa1,qa2],ignore_index=True);cov=pd.concat([cov1,cov2],ignore_index=True)
    wt(qa,'ADDITIVITY_QA.tsv');wt(cov,'CONTRIBUTION_COVERAGE.tsv')
    (OUT/'CONTRIBUTION_DEFINITION.md').write_text('# Additive provenance contribution definition\n\n`C_k=(1/95) × Σ z_g` over frozen members of class k, and `MES=C_MALIGNANT+C_ECOLOGICAL+C_SHARED+C_UNSTABLE`. Classes are never re-standardized or averaged by class size and are not new signatures. Platform-absent genes are structural zero contributions on the within-layer standardized mean-reference scale; coverage is reported. Signed or absolute relative composition summaries are descriptive and are not biological source percentages.\n')
    if (qa.status!='PASS').any(): raise SystemExit('ADDITIVITY_OR_PSEUDOBULK_REPRODUCTION_QA_FAIL')
    (RUN/'session_info/MODULE_A_BULK_PSEUDOBULK_DISCORDANCE_SESSION_INFO.txt').write_text(f'python={platform.python_version()}\nnumpy={np.__version__}\npandas={pd.__version__}\nscipy={__import__("scipy").__version__}\n')
    print('MODULE_A_BULK_PSEUDOBULK_DISCORDANCE_COMPLETE',len(bulk),len(pd.read_csv(OUT/'PSEUDOBULK_FROZEN_CLASS_RESULTS.tsv',sep='\t')))


if __name__=='__main__': main()
