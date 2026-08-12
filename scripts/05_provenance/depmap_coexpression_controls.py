#!/usr/bin/env python3
from __future__ import annotations

import json
import re
from collections import defaultdict
import os
from pathlib import Path

DATA_ROOT = Path(os.environ.get("GBM_MES_DATA_ROOT", "data"))
RESULT_ROOT = Path(os.environ.get("GBM_MES_OUTPUT_ROOT", "results"))
REPO_ROOT = Path(os.environ.get("GBM_MES_REPO_ROOT", Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = RESULT_ROOT / "depmap_null"
V2 = RESULT_ROOT / "depmap"
DATA = DATA_ROOT / "depmap"
SEED = 20260813
TARGET_EFFECTIVE = 500
MAX_ATTEMPTS = 5000
TOL = .20
MIN_PAIRS = 30
NBOOT = 2000
NPERM = 5000
GENE_RE = re.compile(r"\s+\(\d+\)$")


def clean_gene(x: str) -> str:
    return GENE_RE.sub("", str(x)).upper()


def qbins(values: pd.Series, q: int) -> pd.Series:
    return pd.qcut(values.rank(method="first"), q=q, labels=False, duplicates="drop").astype(int)


def greedy_caliper(a: np.ndarray, b: np.ndarray, caliper: float) -> list[tuple[int, int]]:
    cand = sorted((abs(a[i]-b[j]), i, j) for i in range(len(a)) for j in range(len(b)) if abs(a[i]-b[j]) <= caliper)
    ua, ub, out = set(), set(), []
    for _, i, j in cand:
        if i not in ua and j not in ub:
            out.append((i, j)); ua.add(i); ub.add(j)
    return out


def match_ids(ng: list[str], tr: list[str], score: pd.Series) -> tuple[list[tuple[str, str]], float]:
    a, b = score.loc[ng].to_numpy(float), score.loc[tr].to_numpy(float)
    cal = .2 * np.std(np.r_[a, b], ddof=1)
    return [(ng[i], tr[j]) for i, j in greedy_caliper(a, b, cal)], float(cal)


def pair_corr(pairs: list[tuple[str, str]], z: pd.DataFrame, genes: list[str]) -> tuple[float, float, float, float]:
    cors, cosines, ranksims, distances = [], [], [], []
    p = len(genes)
    for a, b in pairs:
        x, y = z.loc[a, genes].to_numpy(float), z.loc[b, genes].to_numpy(float)
        cors.append(spearmanr(x, y).statistic)
        cosines.append(np.dot(x, y) / (np.linalg.norm(x)*np.linalg.norm(y)))
        ranksims.append(1 - np.mean(np.abs(rankdata(x)-rankdata(y)))/(p-1))
        distances.append(np.sqrt(np.mean(((x-y)/p)**2)))
    return float(np.mean(cors)), float(np.mean(cosines)), float(np.mean(ranksims)), float(np.mean(distances))


def rank_norm(x: np.ndarray) -> np.ndarray:
    r = rankdata(x, axis=0, method="average")
    r -= r.mean(axis=0)
    den = np.sqrt(np.sum(r*r, axis=0)); valid = den > 0
    r[:, valid] /= den[valid]; r[:, ~valid] = np.nan
    return r.astype(np.float32)


def summary(norm: np.ndarray, idx: list[int]) -> dict[str, float]:
    c = norm[:, idx].T @ norm[:, idx]
    v = np.abs(c[np.triu_indices(len(idx), 1)]).astype(float)
    q25, med, q75 = np.quantile(v, [.25,.5,.75])
    return {"mean_abs_rho":v.mean(),"median_abs_rho":med,"q25_abs_rho":q25,"q75_abs_rho":q75,"IQR_abs_rho":q75-q25,"fraction_abs_rho_gt_0.2":np.mean(v>.2),"fraction_abs_rho_gt_0.4":np.mean(v>.4)}


def bin_keys(m: int, v: int, d: int, level: int) -> list[tuple[int,int,int]]:
    if level == 0: return [(m,v,d)]
    if level == 1: return [(m,v,dd) for dd in (d-1,d+1) if 0<=dd<=4]
    if level == 2: return [(m,vv,d) for vv in (v-1,v+1) if 0<=vv<=9]
    return [(mm,vv,d) for mm in (m-1,m+1) for vv in (v-1,v+1) if 0<=mm<=9 and 0<=vv<=9]


def construct(template: list[str], features: pd.DataFrame, pools: dict, norm: np.ndarray, gene_idx: dict[str,int], targets: tuple[float,float,float], rng: np.random.Generator) -> tuple[list[str], list[str]]:
    order = rng.permutation(template).tolist()
    selected: list[str] = []
    levels_used: list[str] = []
    used: set[str] = set()
    t25, t50, t75 = targets
    for gene in order:
        row = features.loc[gene]; m,v,d = int(row.mean_bin),int(row.var_bin),int(row.detect_bin)
        candidates, level_name = [], "FAILED"
        for level, name in enumerate(["EXACT_BIN","L1_ADJACENT_DETECTION","L2_ADJACENT_VARIANCE","L3_ADJACENT_MEAN_VARIANCE"]):
            candidates = [g for key in bin_keys(m,v,d,level) for g in pools.get(key,[]) if g not in used]
            candidates = list(dict.fromkeys(candidates))
            if candidates: level_name=name; break
        if not candidates: return [], []
        if len(selected) < 4:
            chosen = str(rng.choice(candidates))
        else:
            ci = [gene_idx[g] for g in candidates]
            si = [gene_idx[g] for g in selected]
            vals = np.abs(norm[:, ci].T @ norm[:, si]).astype(float)
            q25 = np.quantile(vals,.25,axis=1); med=np.quantile(vals,.5,axis=1); q75=np.quantile(vals,.75,axis=1)
            obj = ((q25-t25)/t25)**2 + ((med-t50)/t50)**2 + ((q75-t75)/t75)**2
            top = np.argsort(obj)[:min(5,len(obj))]
            chosen = candidates[int(rng.choice(top))]
        selected.append(chosen); used.add(chosen); levels_used.append(level_name)
    return selected, levels_used


def empirical_lower(null: np.ndarray, obs: float) -> float:
    return float((1+np.sum(null<=obs))/(len(null)+1))


def main() -> None:
    rng=np.random.default_rng(SEED)
    sig=pd.read_csv(V2/'02_signature_reference/signature_registry.tsv',sep='\t')
    scores=pd.read_csv(V2/'03_expression_scores/model_level_scores.tsv',sep='\t').set_index('ModelID')
    ng=scores.index[scores.NextGen_or_traditional.eq('NextGen')].tolist(); tr=scores.index[scores.NextGen_or_traditional.eq('traditional')].tolist(); ids=ng+tr
    hdr=pd.read_csv(DATA/'next_gen_expression.csv',nrows=0).columns.tolist(); dt={c:np.float32 for c in hdr[1:]};dt[hdr[0]]=str
    a=pd.read_csv(DATA/'next_gen_expression.csv',index_col=0,dtype=dt);b=pd.read_csv(DATA/'traditional_expression.csv',index_col=0,dtype=dt)
    raw=pd.concat([a,b]);raw.index=raw.index.astype(str)
    mp={}
    for c in raw.columns: mp.setdefault(clean_gene(c),c)
    cc=list(mp.values());u=raw[cc].copy();u.columns=[clean_gene(c) for c in cc]
    mu=u.mean();sd=u.std(ddof=1); valid=sd[(sd>0)&np.isfinite(sd)].index
    u=u[valid];mu=mu[valid];sd=sd[valid];z=((u.loc[ids]-mu)/sd).astype(np.float32)
    auth=sig.loc[sig.signature_id.eq('CANONICAL_MES95'),'gene'].astype(str).str.upper().drop_duplicates().tolist();mes=[g for g in auth if g in z.columns];assert len(mes)==91
    f=pd.DataFrame({'mean':mu,'variance':u.var(ddof=1),'detection':(u>0).mean()}).dropna();f['mean_bin']=qbins(f['mean'],10);f['var_bin']=qbins(f['variance'],10);f['detect_bin']=qbins(f['detection'],5)
    pools=defaultdict(list)
    for g,r in f.loc[~f.index.isin(set(auth))].iterrows():pools[(int(r.mean_bin),int(r.var_bin),int(r.detect_bin))].append(g)
    ngz=z.loc[ng];trz=z.loc[tr];res=pd.concat([(ngz-ngz.mean())/ngz.std(ddof=1),(trz-trz.mean())/trz.std(ddof=1)])
    norm=rank_norm(res.to_numpy(float));gidx={g:i for i,g in enumerate(z.columns)}
    target=summary(norm,[gidx[g] for g in mes]);targets=(target['q25_abs_rho'],target['median_abs_rho'],target['q75_abs_rho'])
    mes_obs=pd.read_csv(ROOT/'02_exact91_random_null/mes95_exact91_observed.tsv',sep='\t').iloc[0];mes_corr=float(mes_obs.contribution_correlation)
    accepted=[];check=[]
    for attempt in range(1,MAX_ATTEMPTS+1):
        genes,levels=construct(mes,f,pools,norm,gidx,targets,rng)
        if not genes:
            check.append({'attempt':attempt,'construction_success':False,'coexpression_match':False,'pair_eligible':False});continue
        cx=summary(norm,[gidx[g] for g in genes])
        coex_ok=abs(cx['median_abs_rho']/target['median_abs_rho']-1)<=TOL and abs(cx['IQR_abs_rho']/target['IQR_abs_rho']-1)<=TOL
        if not coex_ok:
            check.append({'attempt':attempt,'construction_success':True,'coexpression_match':False,'pair_eligible':False,**cx});continue
        score=z[genes].mean(axis=1);pairs,cal=match_ids(ng,tr,score);pair_ok=len(pairs)>=MIN_PAIRS
        check.append({'attempt':attempt,'construction_success':True,'coexpression_match':True,'pair_eligible':pair_ok,'n_matched_pairs':len(pairs),**cx})
        if not pair_ok: continue
        pc,cos,rs,dist=pair_corr(pairs,z,genes)
        accepted.append({'signature_id':len(accepted)+1,'attempt':attempt,'n_genes':91,'n_exact_bin':sum(x=='EXACT_BIN' for x in levels),'n_near_bin':sum(x!='EXACT_BIN' for x in levels),'n_matched_pairs':len(pairs),'caliper':cal,'contribution_correlation':pc,'cosine_similarity':cos,'rank_similarity':rs,'contribution_distance':dist,**{f'coexpression_{k}':v for k,v in cx.items()},'genes':';'.join(genes),'matching_levels':';'.join(levels)})
        if len(accepted)>=TARGET_EFFECTIVE:break
    acc=pd.DataFrame(accepted);pd.DataFrame(check).to_csv(ROOT/'03_coexpression_null/coexpression_greedy_attempt_check.tsv.gz',sep='\t',index=False,compression='gzip')
    acc.to_csv(ROOT/'03_coexpression_null/coexpression_null.tsv',sep='\t',index=False)
    vals=acc.contribution_correlation.to_numpy(float) if len(acc) else np.array([])
    p=empirical_lower(vals,mes_corr) if len(vals) else np.nan; pct=100*np.mean(vals>=mes_corr) if len(vals) else np.nan
    decision='COEXPRESSION_NULL_LOW_POWER' if len(acc)<500 else ('COEXPRESSION_AWARE_NULL_PASS' if pct>=95 and p<.05 else 'COEXPRESSION_AWARE_NULL_FAIL')
    cs=pd.DataFrame([{'construction_method':'outcome-blind greedy q25/median/q75 matching within expression bins','attempts_used':int(pd.DataFrame(check).attempt.max()),'effective_random_signatures_n':len(acc),'chosen_relative_tolerance':TOL,'MES95_median_abs_rho':target['median_abs_rho'],'MES95_IQR_abs_rho':target['IQR_abs_rho'],'MES95_contribution_correlation':mes_corr,'coexpression_null_mean':vals.mean() if len(vals) else np.nan,'coexpression_null_SD':vals.std(ddof=1) if len(vals)>1 else np.nan,'MES95_divergence_percentile':pct,'empirical_P':p,'decision':decision,'selection_did_not_use_outcome':True}])
    cs.to_csv(ROOT/'03_coexpression_null/coexpression_summary.tsv',sep='\t',index=False)

    # Complete the three frozen supplementary controls interrupted in phase 1.
    controls=[]
    for sid in ['NEFTEL_AC_LIKE','HALLMARK_EMT','NEFTEL_NPC_LIKE']:
        ag=sig.loc[sig.signature_id.eq(sid),'gene'].astype(str).str.upper().drop_duplicates().tolist();genes=[g for g in ag if g in z.columns]
        score=z[genes].mean(axis=1);pairs,cal=match_ids(ng,tr,score);obs=pair_corr(pairs,z,genes)[0]
        pairvals=[]
        for x,y in pairs:pairvals.append(spearmanr(z.loc[x,genes],z.loc[y,genes]).statistic)
        boots=np.empty(NBOOT)
        for i in range(NBOOT):boots[i]=np.mean(rng.choice(pairvals,len(pairvals),replace=True))
        perm=np.empty(NPERM);arr=np.array(ids)
        for i in range(NPERM):
            sh=rng.permutation(arr);pp,_=match_ids(sh[:len(ng)].tolist(),sh[len(ng):].tolist(),score);perm[i]=pair_corr(pp,z,genes)[0]
        controls.append({'signature_id':sid,'reference_gene_n':len(ag),'measurable_gene_n':len(genes),'matched_pair_n':len(pairs),'caliper':cal,'contribution_correlation':obs,'bootstrap_CI_low':np.quantile(boots,.025),'bootstrap_CI_high':np.quantile(boots,.975),'permutation_null_mean':perm.mean(),'empirical_P_lower_similarity':empirical_lower(perm,obs),'role':'SUPPLEMENTARY_CONTROL'})
    ctrl=pd.DataFrame(controls);ctrl.to_csv(ROOT/'04_cross_null_comparison/positive_control_signatures.tsv',sep='\t',index=False)
    cs.to_csv(ROOT/'06_tables/coexpression_summary.tsv',sep='\t',index=False);acc.drop(columns=['genes','matching_levels'],errors='ignore').to_csv(ROOT/'06_tables/coexpression_null.tsv',sep='\t',index=False);ctrl.to_csv(ROOT/'06_tables/positive_control_signatures.tsv',sep='\t',index=False)
    receipt={'seed':SEED,'construction_target':TARGET_EFFECTIVE,'max_attempts':MAX_ATTEMPTS,'effective_n':len(acc),'tolerance':TOL,'decision':decision,'outcome_used_in_construction':False,'dependency_branch':'NOT_TOUCHED_FROZEN'}
    (ROOT/'10_logs/phase2_coexpression_greedy_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':main()
