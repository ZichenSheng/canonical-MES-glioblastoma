#!/usr/bin/env python3
from __future__ import annotations
import importlib.util, json, math, platform
import os
from pathlib import Path

DATA_ROOT = Path(os.environ.get("GBM_MES_DATA_ROOT", "data"))
RESULT_ROOT = Path(os.environ.get("GBM_MES_OUTPUT_ROOT", "results"))
REPO_ROOT = Path(os.environ.get("GBM_MES_REPO_ROOT", Path(__file__).resolve().parents[2]))
import numpy as np, pandas as pd
from scipy import stats
from scipy.optimize import minimize_scalar

RUN=RESULT_ROOT/'provenance';D=RUN/'ablation';A=RUN/'contributions';V13=RESULT_ROOT/'interpretability';V11=RESULT_ROOT/'scoring';V1=DATA_ROOT/'prepared/core';SIG=REPO_ROOT/'resources/signatures/signature_gene_sets.tsv'
CLASSES=['MALIGNANT_DOMINANT','ECOLOGICAL_DOMINANT','SHARED_MIXED_ORIGIN','LOW_INFORMATION_OR_UNSTABLE'];SEED=20260801406

def wt(x,n): x.to_csv(D/n,sep='\t',index=False,na_rep='NA')
def spcor_cols(Y,x):
    ry=stats.rankdata(Y,axis=0);rx=stats.rankdata(x);ry-=ry.mean(axis=0);rx=rx-rx.mean();return (ry*rx[:,None]).sum(axis=0)/np.sqrt((ry*ry).sum(axis=0)*(rx*rx).sum())
def reml(yi,vi):
    yi=np.asarray(yi,float);vi=np.asarray(vi,float)
    def obj(t):w=1/(vi+t);mu=np.sum(w*yi)/sum(w);return .5*(np.sum(np.log(vi+t))+np.log(sum(w))+np.sum(w*(yi-mu)**2))
    tau=max(0,float(minimize_scalar(obj,bounds=(0,10),method='bounded').x));w=1/(vi+tau);mu=float(np.sum(w*yi)/sum(w));se=math.sqrt(1/sum(w));return mu,mu-1.96*se,mu+1.96*se,mu-1.96*math.sqrt(tau+se*se),mu+1.96*math.sqrt(tau+se*se),tau
def huber_residuals(Y,m,it=15):
    m=np.asarray(m,float);X=np.column_stack([np.ones(len(m)),m]);b=np.linalg.pinv(X)@Y
    for _ in range(it):
        r=Y-X@b;scale=1.4826*np.median(np.abs(r-np.median(r,axis=0)),axis=0)+1e-8;w=np.minimum(1,1.345*scale[None,:]/(np.abs(r)+1e-12));sw=w.sum(axis=0);sm=(w*m[:,None]).sum(axis=0);smm=(w*(m*m)[:,None]).sum(axis=0);sy=(w*Y).sum(axis=0);smy=(w*m[:,None]*Y).sum(axis=0);den=sw*smm-sm*sm+1e-12;b1=(sw*smy-sm*sy)/den;b0=(sy-b1*sm)/sw;b=np.vstack([b0,b1])
    return Y-X@b

def load_freeze():
    p=pd.read_csv(RUN/'01_provenance_freeze/PROVENANCE_CLASS_FREEZE.tsv',sep='\t');g='gene' if 'gene'in p else 'harmonized_symbol';c='primary_origin' if 'primary_origin'in p else 'origin_class';p[g]=p[g].str.upper();return p[g].tolist(),dict(zip(p[g],p[c]))

def raw_feature_table(genes,cmap):
    spec=importlib.util.spec_from_file_location('ma',Path(__file__).with_name('provenance_contributions.py'));ma=importlib.util.module_from_spec(spec);spec.loader.exec_module(ma)
    feats=[]
    for co,path,mfile,lg in [('CGGA325',DATA_ROOT/'bulk/cgga325_expression.tsv','membership_cgga325.tsv',True),('CGGA693',DATA_ROOT/'bulk/cgga693_expression.tsv','membership_cgga693.tsv',True),('TCGA',DATA_ROOT/'bulk/tcga_gbm_expression.tsv.gz','membership_tcga.tsv',False)]:
        md=pd.read_csv(DATA_ROOT/'prepared/cohorts'/mfile,sep='\t');ids=md.entity_id.astype(str).tolist()
        if co=='TCGA':
            tm=pd.read_csv(V11/'05_tcga_strict/TCGA_STRICT_COHORT.tsv',sep='\t');ids=tm[tm.strict_inclusion.astype(bool)&tm.patient_id.astype(str).isin(ids)].sample_id.astype(str).tolist()
        x=ma.read_bulk_matrix(path,ids,lg).reindex(genes);mu=x.mean(axis=1);vr=np.log1p(x.var(axis=1));mu=(mu-mu.mean())/mu.std();vr=(vr-vr.mean())/vr.std();feats.append(pd.DataFrame({'gene':genes,f'mean_{co}':mu.values,f'var_{co}':vr.values}))
    f=feats[0]
    for q in feats[1:]:f=f.merge(q,on='gene')
    f['bulk_mean_expression']=f[[c for c in f if c.startswith('mean_')]].mean(axis=1);f['bulk_expression_variance']=f[[c for c in f if c.startswith('var_')]].mean(axis=1)
    con=pd.read_csv(V13/'03_gene_provenance/GENE_SCRNA_SNRNA_CONCORDANCE.tsv',sep='\t');con=con[con.metric.eq('detection_fraction')][['gene','GBmap','GSE174554']].rename(columns={'GBmap':'scRNA_detection','GSE174554':'snRNA_detection'})
    cov=pd.read_csv(V13/'03_gene_provenance/CANONICAL_95_COVERAGE.tsv',sep='\t').rename(columns={'harmonized_symbol':'gene'});pm=pd.read_csv(V13/'03_gene_provenance/CANONICAL_95_PROVENANCE_MAP.tsv',sep='\t')[['gene','malignant_median_contribution']]
    f=f[['gene','bulk_mean_expression','bulk_expression_variance']].merge(con,on='gene').merge(cov[['gene','GBmap_present','GSE174554_present']],on='gene').merge(pm,on='gene')
    f['platform_coverage']=f.GBmap_present.astype(int)+f.GSE174554_present.astype(int)+f.malignant_median_contribution.notna().astype(int)
    f['scRNA_detection_missing']=f.scRNA_detection.isna().astype(int);f['snRNA_detection_missing']=f.snRNA_detection.isna().astype(int)
    f['scRNA_detection']=f.scRNA_detection.fillna(f.scRNA_detection.median());f['snRNA_detection']=f.snRNA_detection.fillna(f.snRNA_detection.median())
    reg=pd.read_csv(SIG,sep='\t');reg['gene']=reg.gene_symbol_clean.astype(str).str.strip().str.upper();ctx=set(reg.loc[reg.signature_id.isin(['HALLMARK_HYPOXIA','HALLMARK_EMT','NABA_CORE_MATRISOME','MYELOID_CORE_IDENTITY']),'gene']);f['direct_overlap_status']=f.gene.isin(ctx).astype(int)
    md=f[['scRNA_detection','snRNA_detection']].min(axis=1);f['dropout_low']=(md<.1).astype(int);f['dropout_mid']=((md>=.1)&(md<.3)).astype(int);f['dropout_high']=(md>=.3).astype(int);f['origin_class']=f.gene.map(cmap);return f

def matched_sets(features,genes,cmap):
    cols=['bulk_mean_expression','bulk_expression_variance','scRNA_detection','snRNA_detection','scRNA_detection_missing','snRNA_detection_missing','platform_coverage','direct_overlap_status','dropout_low','dropout_mid','dropout_high'];X=features.set_index('gene').loc[genes,cols].to_numpy(float);X=(X-X.mean(axis=0))/(X.std(axis=0,ddof=1)+1e-12);rng=np.random.default_rng(SEED);allrows=[];diag=[]
    for ci,cl in enumerate(CLASSES):
        target=np.array([i for i,g in enumerate(genes) if cmap[g]==cl]);k=len(target);tm=X[target].mean(axis=0);seen={};attempts=0
        while len(seen)<120000 and attempts<400000:
            idx=tuple(sorted(rng.choice(len(genes),k,replace=False).tolist()));attempts+=1
            if idx==tuple(target.tolist()) or idx in seen:continue
            diff=X[list(idx)].mean(axis=0)-tm;mx=float(np.max(np.abs(diff)));dist=float(np.sqrt(np.mean(diff**2)));level='STRICT_0.25' if mx<=.25 else 'RELAXED_0.35' if mx<=.35 else 'RELAXED_0.50' if mx<=.50 else 'RANKED_FALLBACK';seen[idx]=(mx,dist,level,diff)
        cand=sorted(seen.items(),key=lambda kv:({'STRICT_0.25':0,'RELAXED_0.35':1,'RELAXED_0.50':2,'RANKED_FALLBACK':3}[kv[1][2]],kv[1][1]))[:5000]
        if len(cand)<5000:raise RuntimeError(f'insufficient matched sets {cl}')
        for rid,(idx,(mx,dist,lev,diff)) in enumerate(cand,1):allrows.append(dict(target_class=cl,null_replicate=rid,deleted_gene_count=k,deleted_genes=';'.join(genes[i] for i in idx),matching_level=lev,max_abs_standardized_imbalance=mx,matching_distance=dist))
        arr=np.array([list(x[0]) for x in cand]);ach=X[arr].mean(axis=1).mean(axis=0)
        for j,col in enumerate(cols):diag.append(dict(target_class=cl,feature=col,target_mean_standardized=tm[j],null_mean_standardized=ach[j],standardized_difference=ach[j]-tm[j],n_null=5000,strict_n=sum(x[1][2]=='STRICT_0.25' for x in cand),relaxed_035_n=sum(x[1][2]=='RELAXED_0.35' for x in cand),relaxed_050_n=sum(x[1][2]=='RELAXED_0.50' for x in cand),fallback_n=sum(x[1][2]=='RANKED_FALLBACK' for x in cand),candidate_attempts=attempts))
    sets=pd.DataFrame(allrows);wt(sets,'MATCHED_DELETION_SETS_INTERNAL.tsv');wt(pd.DataFrame(diag),'MATCHING_DIAGNOSTICS.tsv');wt(features,'MATCHING_FEATURES_INTERNAL.tsv');return sets

def deletion_matrix(rows,genes):
    M=np.zeros((len(rows),len(genes)),float);gm={g:i for i,g in enumerate(genes)}
    for i,s in enumerate(rows.deleted_genes):M[i,[gm[g] for g in s.split(';')]]=1
    return M

def pseudo_weights(design):
    out={}
    for (exp,comp),ii0 in design.groupby(['experiment','compartment']).groups.items():
        if (exp=='A_FIXED_MALIGNANT' and comp not in ['myeloid','endothelial','pericyte']) or (exp=='C_TWO_DIMENSIONAL' and comp not in ['myeloid','endothelial','pericyte']):continue
        q=design.loc[ii0];pat=pd.get_dummies(q.patient_key,drop_first=True,dtype=float).to_numpy();mes=q.malignant_mes_fraction.to_numpy(float);com=q.composition_fraction.to_numpy(float)
        if exp=='A_FIXED_MALIGNANT':terms=np.column_stack([com]);names=['composition_slope']
        elif exp=='B_FIXED_COMPOSITION':terms=np.column_stack([mes]);names=['malignant_state_slope']
        else:terms=np.column_stack([mes,com,mes*com]);names=['malignant_state_slope','composition_slope','interaction']
        X=np.column_stack([np.ones(len(q)),terms,pat]);P=np.linalg.pinv(X)
        for j,nm in enumerate(names,1):out[f'PSEUDOBULK_{exp}_{comp}_{nm.upper()}']=(np.asarray(list(ii0)),P[j])
    return out

def strict_endpoints(Zs,ids,context,deleted,genes,overlap):
    vals={};meta_raw=[];meta_pr=[];ns=[]
    for co,Z in Zs.items():
        total=Z.sum(axis=1)/95;abl=total-Z[:,deleted].sum(axis=1)/95;ctx=context.set_index(['cohort','patient_id']).loc[[(co,x) for x in ids[co]],'context'].to_numpy(float);r0=stats.spearmanr(total,ctx).statistic;r1=stats.spearmanr(abl,ctx).statistic;vals[f'STRICT_CONTEXT_RHO_{co}_RAW']=(r0,r1,None)
        keep=[j for j,g in enumerate(genes) if g not in overlap];tpr=Z[:,keep].sum(axis=1)/95;apr=tpr-Z[:,[j for j in deleted if j in keep]].sum(axis=1)/95;p0=stats.spearmanr(tpr,ctx).statistic;p1=stats.spearmanr(apr,ctx).statistic;vals[f'STRICT_CONTEXT_RHO_{co}_OVERLAP_PRUNED']=(p0,p1,None);meta_raw.append((r0,r1));meta_pr.append((p0,p1));ns.append(len(ctx))
    for label,arr in [('RAW',meta_raw),('OVERLAP_PRUNED',meta_pr)]:
        vi=1/(np.asarray(ns)-3);m0=reml(np.arctanh(np.clip([x[0] for x in arr],-.999999,.999999)),vi);m1=reml(np.arctanh(np.clip([x[1] for x in arr],-.999999,.999999)),vi);vals[f'STRICT_CONTEXT_META_RHO_{label}']=(np.tanh(m0[0]),np.tanh(m1[0]),(np.tanh(m0[3]),np.tanh(m0[4]),np.tanh(m1[3]),np.tanh(m1[4])))
    return vals

def discord_endpoints(Z,ids,dis,deleted_sets,co):
    q=dis[dis.cohort.eq(co)].set_index('patient_id').loc[ids];total=Z.sum(axis=1)/95;Y=total[:,None]-Z@deleted_sets.T/95;Y=(Y-Y.mean(axis=0))/(Y.std(axis=0,ddof=1)+1e-12);m=q.M.to_numpy(float);ctx=q['context'].to_numpy(float);rho=spcor_cols(Y,m);qtM=m>=np.quantile(m,.75);thr=np.quantile(Y,.75,axis=0);false=((Y>=thr[None,:])&(~qtM[:,None])).mean(axis=0);res=huber_residuals(Y,m);rsd=res.std(axis=0,ddof=1);eco=spcor_cols(res,ctx);return {f'DISCORDANCE_BULK_MALIGNANT_RHO_{co}':rho,f'DISCORDANCE_TOPQ_FALSE_HIGH_{co}':false,f'DISCORDANCE_RESIDUAL_SD_{co}':rsd,f'DISCORDANCE_ECOLOGY_RHO_{co}':eco}

def main():
    genes,cmap=load_freeze();features=raw_feature_table(genes,cmap);sets=matched_sets(features,genes,cmap)
    npz=np.load(A/'STRICT_BULK_GENE_Z_INTERNAL.npz',allow_pickle=True);Zs={co:npz[f'Z_{co}'] for co in ['CGGA325','CGGA693','TCGA']};ids={co:npz[f'IDS_{co}'].astype(str).tolist() for co in Zs};context=pd.read_csv(A/'STRICT_BULK_CLASS_CONTRIBUTIONS.tsv',sep='\t');reg=pd.read_csv(SIG,sep='\t');reg['gene']=reg.gene_symbol_clean.astype(str).str.strip().str.upper();overlap=set(reg.loc[reg.signature_id.isin(['HALLMARK_HYPOXIA','HALLMARK_EMT','NABA_CORE_MATRISOME','MYELOID_CORE_IDENTITY']),'gene'])
    pn=np.load(A/'PSEUDOBULK_GENE_Z_INTERNAL.npz',allow_pickle=True);PZ=pn['Z'];design=pd.read_csv(V11/'04_pseudobulk_benchmark/PSEUDOBULK_DESIGN.tsv',sep='\t').set_index('pseudobulk_id').loc[pn['pseudobulk_id'].astype(str)].reset_index();pweights=pseudo_weights(design);dis=pd.read_csv(V13/'02_patient_discordance/PATIENT_LEVEL_DISCORDANCE.tsv',sep='\t')
    results=[];rescaled=[];nullrows=[]
    for cl in CLASSES:
        true=np.array([i for i,g in enumerate(genes) if cmap[g]==cl]);qsets=sets[sets.target_class.eq(cl)].sort_values('null_replicate');Mnull=deletion_matrix(qsets,genes);Mall=np.vstack([np.zeros(95),np.eye(95)[true].sum(axis=0),Mnull])
        # Strict endpoints are nonlinear and evaluated set-by-set.
        originals=strict_endpoints(Zs,ids,context,[],genes,overlap)
        truevals=strict_endpoints(Zs,ids,context,true.tolist(),genes,overlap)
        for ep,(o,a,aux) in truevals.items():
            oo=originals[ep][0];row=dict(target_class=cl,endpoint=ep,original_value=oo,ablated_value=a,ablation_effect=a-oo,deleted_gene_count=len(true),primary_denominator=95,interpretation='measurement behavior contribution; not causal source')
            if aux is not None:row.update(original_prediction_low=aux[0],original_prediction_high=aux[1],ablated_prediction_low=aux[2],ablated_prediction_high=aux[3])
            results.append(row)
        for rid,deleted in enumerate(np.where(Mnull>0)[1].reshape(len(Mnull),-1),1):
            vv=strict_endpoints(Zs,ids,context,deleted.tolist(),genes,overlap)
            for ep,(o,a,aux) in vv.items():nullrows.append(dict(target_class=cl,null_replicate=rid,endpoint=ep,original_value=originals[ep][0],ablated_value=a,ablation_effect=a-originals[ep][0],matching_level=qsets.iloc[rid-1].matching_level))
        # Pseudo-bulk slopes are linear in the score.
        total=PZ.sum(axis=1)/95
        for ep,(ii,w) in pweights.items():
            gene_eff=(w[:,None]*PZ[ii]).sum(axis=0)/95;o=float(w@total[ii]);a=o-gene_eff[true].sum();results.append(dict(target_class=cl,endpoint=ep,original_value=o,ablated_value=a,ablation_effect=a-o,deleted_gene_count=len(true),primary_denominator=95,interpretation='measurement behavior contribution; not causal source'))
            ar=(o-gene_eff[true].sum())*95/(95-len(true));rescaled.append(dict(target_class=cl,endpoint=ep,original_value=o,ablated_value=ar,ablation_effect=ar-o,remaining_gene_denominator=95-len(true)))
            eff=-(Mnull@gene_eff)
            for rid,e in enumerate(eff,1):nullrows.append(dict(target_class=cl,null_replicate=rid,endpoint=ep,original_value=o,ablated_value=o+e,ablation_effect=e,matching_level=qsets.iloc[rid-1].matching_level))
        # Discordance endpoints: column 0 original, 1 true, then nulls.
        for co in ['CGGA325','CGGA693']:
            dd=discord_endpoints(Zs[co],ids[co],dis,Mall,co)
            for ep,v in dd.items():
                results.append(dict(target_class=cl,endpoint=ep,original_value=v[0],ablated_value=v[1],ablation_effect=v[1]-v[0],deleted_gene_count=len(true),primary_denominator=95,interpretation='recomputed robust residual behavior; not causal source'))
                for rid,e in enumerate(v[2:]-v[0],1):nullrows.append(dict(target_class=cl,null_replicate=rid,endpoint=ep,original_value=v[0],ablated_value=v[rid+1],ablation_effect=e,matching_level=qsets.iloc[rid-1].matching_level))
        # Primary/rescaled strict correlations are scale invariant; record true rows.
        for ep,(o,a,aux) in truevals.items():rescaled.append(dict(target_class=cl,endpoint=ep,original_value=o,ablated_value=a,ablation_effect=a-o,remaining_gene_denominator=95-len(true)))
    wt(pd.DataFrame(results),'CLASS_ABLATION_RESULTS.tsv');wt(pd.DataFrame(rescaled),'CLASS_ABLATION_RESCALED_SENSITIVITY.tsv');wt(pd.DataFrame(nullrows),'MATCHED_DELETION_NULL.tsv')
    (RUN/'session_info/MODULE_B_ABLATION_CORE_SESSION_INFO.txt').write_text(f'python={platform.python_version()}\nnumpy={np.__version__}\npandas={pd.__version__}\nscipy={__import__("scipy").__version__}\n')
    print('MODULE_B_ABLATION_CORE_COMPLETE',len(results),len(nullrows),len(sets))
if __name__=='__main__':main()
