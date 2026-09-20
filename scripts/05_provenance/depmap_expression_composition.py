#!/usr/bin/env python3
"""Frozen signature scoring, context effects, MES95 contributions, and same-score matching."""

from __future__ import annotations

import json
import math
import re
import os
from pathlib import Path

DATA_ROOT = Path(os.environ.get("GBM_MES_DATA_ROOT", "data"))
RESULT_ROOT = Path(os.environ.get("GBM_MES_OUTPUT_ROOT", "results"))
REPO_ROOT = Path(os.environ.get("GBM_MES_REPO_ROOT", Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cosine
from scipy.stats import mannwhitneyu, pearsonr, rankdata, spearmanr


ROOT = RESULT_ROOT / "depmap"
DATA = DATA_ROOT / "depmap"
PATIENT_REF = RESULT_ROOT / "interpretability" / "03_gene_provenance" / "CANONICAL_95_PROVENANCE_MAP.tsv"
SEED = 20260809
NBOOT = 2000
NPERM = 2000
GENE_RE = re.compile(r"\s+\(\d+\)$")


def clean_gene(x: str) -> str:
    return GENE_RE.sub("", str(x)).upper()


def hedges_g(a: np.ndarray, b: np.ndarray) -> float:
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    sp = math.sqrt(((len(a)-1)*np.var(a, ddof=1) + (len(b)-1)*np.var(b, ddof=1)) / (len(a)+len(b)-2))
    if sp == 0:
        return 0.0
    d = (np.mean(a)-np.mean(b))/sp
    return d * (1 - 3/(4*(len(a)+len(b))-9))


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    return float((np.greater(a[:, None], b).sum() - np.less(a[:, None], b).sum()) / (len(a)*len(b)))


def boot_diff(a: np.ndarray, b: np.ndarray, rng: np.random.Generator, fun=np.median, nboot=NBOOT) -> tuple[float, float, float]:
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    vals = np.empty(nboot)
    for i in range(nboot):
        vals[i] = fun(rng.choice(a, len(a), replace=True)) - fun(rng.choice(b, len(b), replace=True))
    return float(fun(a)-fun(b)), float(np.quantile(vals, .025)), float(np.quantile(vals, .975))


def corr_safe(a: np.ndarray, b: np.ndarray, kind="spearman") -> float:
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3:
        return np.nan
    return float(spearmanr(a[m], b[m]).statistic if kind == "spearman" else pearsonr(a[m], b[m]).statistic)


def greedy_caliper(a: np.ndarray, b: np.ndarray, caliper: float) -> list[tuple[int, int]]:
    cand = sorted((abs(a[i]-b[j]), i, j) for i in range(len(a)) for j in range(len(b)) if abs(a[i]-b[j]) <= caliper)
    used_a, used_b, out = set(), set(), []
    for _, i, j in cand:
        if i not in used_a and j not in used_b:
            out.append((i, j)); used_a.add(i); used_b.add(j)
    return out


def pair_metrics(ng_ids, tr_ids, contrib: pd.DataFrame, components: pd.DataFrame, scores: pd.Series, method: str) -> pd.DataFrame:
    rows = []
    for k, (ng, tr) in enumerate(zip(ng_ids, tr_ids), 1):
        x, y = contrib.loc[ng].to_numpy(float), contrib.loc[tr].to_numpy(float)
        px, py = components.loc[ng].to_numpy(float), components.loc[tr].to_numpy(float)
        rows.append({
            "method": method, "pair_id": f"{method}_{k:03d}", "NextGen_ModelID": ng, "traditional_ModelID": tr,
            "NextGen_MES95": scores.loc[ng], "traditional_MES95": scores.loc[tr],
            "absolute_score_difference": abs(scores.loc[ng]-scores.loc[tr]),
            "contribution_spearman": corr_safe(x, y, "spearman"),
            "contribution_cosine_similarity": 1-cosine(x, y) if np.linalg.norm(x)>0 and np.linalg.norm(y)>0 else np.nan,
            "contribution_rank_distance": float(np.mean(np.abs(rankdata(x)-rankdata(y)))/(len(x)-1)),
            "provenance_euclidean_distance": float(np.linalg.norm(px-py)),
            **{f"delta_{c}_NextGen_minus_traditional": px[i]-py[i] for i, c in enumerate(components.columns)},
        })
    return pd.DataFrame(rows)


def main() -> None:
    rng = np.random.default_rng(SEED)
    reg = pd.read_csv(ROOT / "01_analysis_units/model_registry.tsv", sep="\t")
    reg = reg.loc[reg["strict_GBM"].eq(True)].set_index("ModelID")
    sig = pd.read_csv(ROOT / "02_signature_reference/signature_registry.tsv", sep="\t")
    pref = pd.read_csv(PATIENT_REF, sep="\t")
    pref["gene"] = pref["gene"].str.upper()

    # Read complete matrices as float32 so within-sample rank scores are computed against the full gene universe.
    exp_cols = pd.read_csv(DATA / "next_gen_expression.csv", nrows=0).columns.tolist()
    exp_types = {c: np.float32 for c in exp_cols[1:]}
    exp_types[exp_cols[0]] = str
    ng = pd.read_csv(DATA / "next_gen_expression.csv", index_col=0, dtype=exp_types)
    tr = pd.read_csv(DATA / "traditional_expression.csv", index_col=0, dtype=exp_types)
    assert list(ng.columns) == list(tr.columns)
    raw = pd.concat([ng, tr], axis=0)
    raw.index = raw.index.astype(str)
    symbols = [clean_gene(c) for c in raw.columns]
    sym_to_cols: dict[str, list[str]] = {}
    for c, s in zip(raw.columns, symbols):
        sym_to_cols.setdefault(s, []).append(c)
    duplicate_symbols = {k: v for k, v in sym_to_cols.items() if len(v) > 1}
    mapped_col = {g: cs[0] for g, cs in sym_to_cols.items()}

    mu = raw.mean(axis=0)
    sd = raw.std(axis=0, ddof=1).replace(0, np.nan)
    z = (raw-mu)/sd
    # Rank percentiles centered to approximately [-0.5, 0.5].
    rr = raw.rank(axis=1, method="average", pct=True)-0.5

    score = pd.DataFrame(index=raw.index)
    rank_score = pd.DataFrame(index=raw.index)
    coverage_rows = []
    for sid, d in sig.groupby("signature_id", sort=False):
        genes = d["gene"].astype(str).str.upper().tolist()
        avail = [g for g in genes if g in mapped_col]
        cols = [mapped_col[g] for g in avail]
        meta = d.drop_duplicates("gene").set_index(d["gene"].astype(str).str.upper())
        weights = np.array([float(meta.loc[g, "weight"])*float(meta.loc[g, "direction"]) for g in avail])
        denom = np.sum(np.abs(weights))
        score[sid] = np.nansum(z[cols].to_numpy()*weights, axis=1)/denom
        rank_score[sid] = np.nansum(rr[cols].to_numpy()*weights, axis=1)/denom
        coverage_rows.append({"signature_id": sid, "reference_gene_n": len(set(genes)), "mapped_gene_n": len(set(avail)),
                              "coverage": len(set(avail))/len(set(genes)), "status": "PASS" if len(set(avail))/len(set(genes)) >= .8 else "LOW_COVERAGE"})
    pd.DataFrame(coverage_rows).to_csv(ROOT / "03_expression_scores/signature_coverage.tsv", sep="\t", index=False)
    pd.DataFrame([{"gene_symbol": k, "columns": ";".join(v), "resolution": "first_column_no_silent_collapse"} for k,v in duplicate_symbols.items()])\
        .to_csv(ROOT / "03_expression_scores/duplicate_gene_symbol_check.tsv", sep="\t", index=False)

    # Model-level score artifact includes only the frozen strict-GBM universe with available expression.
    mids = [x for x in reg.index if x in score.index]
    out = score.loc[mids].copy()
    out.columns = [f"score_z__{c}" for c in out.columns]
    rs = rank_score.loc[mids].copy(); rs.columns = [f"score_rank__{c}" for c in rs.columns]
    out = reg.loc[mids, ["NextGen_or_traditional", "supp57_GBM", "RNA_available", "CRISPR_available", "CN_available", "RNA_CRISPR", "RNA_CRISPR_CN", "Library", "GrowthPattern", "PlateCoating"]].join(out).join(rs)
    out.reset_index().to_csv(ROOT / "03_expression_scores/model_level_scores.tsv", sep="\t", index=False)

    # Exclusive provenance contributions sum to canonical MES95; modality is an overlapping flag and is separate.
    mes = sig.loc[sig["signature_id"].eq("CANONICAL_MES95"), "gene"].str.upper().tolist()
    mes_avail = [g for g in mes if g in mapped_col]
    mes_z = pd.DataFrame({g: z[mapped_col[g]]/len(mes_avail) for g in mes_avail}, index=z.index)
    prov = pref.set_index("gene")["primary_origin"].reindex(mes_avail)
    exclusive = ["MALIGNANT_DOMINANT", "ECOLOGICAL_DOMINANT", "SHARED_MIXED_ORIGIN", "LOW_INFORMATION_OR_UNSTABLE"]
    comp = pd.DataFrame(index=z.index)
    for c in exclusive:
        gg = prov.index[prov.eq(c)].tolist()
        comp[c] = mes_z[gg].sum(axis=1)
    mod_genes = set(sig.loc[sig["signature_id"].eq("MES95_MODALITY_DEPENDENT_PROVENANCE"), "gene"].str.upper())
    comp["MODALITY_DEPENDENT_OVERLAP"] = mes_z[[g for g in mes_avail if g in mod_genes]].sum(axis=1)
    comp.loc[mids].reset_index(names="ModelID").to_csv(ROOT / "03_expression_scores/provenance_component_scores.tsv", sep="\t", index=False)

    # Context distribution effects for all frozen scores and contribution components.
    context_rows = []
    ng_ids = reg.index[(reg["NextGen_or_traditional"].eq("NextGen")) & reg.index.isin(score.index)].tolist()
    tr_ids = reg.index[(reg["NextGen_or_traditional"].eq("traditional")) & reg.index.isin(score.index)].tolist()
    analyses = {f"score_z__{c}": score[c] for c in score.columns}
    analyses.update({f"provenance_contribution__{c}": comp[c] for c in comp.columns})
    for aid, s in analyses.items():
        a, b = s.loc[ng_ids].to_numpy(float), s.loc[tr_ids].to_numpy(float)
        diff, lo, hi = boot_diff(a, b, rng)
        try: p = mannwhitneyu(a, b, alternative="two-sided").pvalue
        except ValueError: p = np.nan
        pooled = np.r_[a,b]
        full_direction = np.sign(np.median(a)-np.median(b))
        loo = []
        for i in range(len(a)): loo.append(np.sign(np.median(np.delete(a,i))-np.median(b)) == full_direction)
        for i in range(len(b)): loo.append(np.sign(np.median(a)-np.median(np.delete(b,i))) == full_direction)
        context_rows.append({"analysis": aid, "NextGen_n": len(a), "traditional_n": len(b),
                             "NextGen_median": np.median(a), "NextGen_IQR": np.quantile(a,.75)-np.quantile(a,.25),
                             "traditional_median": np.median(b), "traditional_IQR": np.quantile(b,.75)-np.quantile(b,.25),
                             "median_difference_NextGen_minus_traditional": diff, "bootstrap_CI_low": lo, "bootstrap_CI_high": hi,
                             "hedges_g": hedges_g(a,b), "cliffs_delta": cliffs_delta(a,b), "MWU_p": p,
                             "LOO_direction_fraction": np.mean(loo), "supported_CI_excludes_zero": bool(lo>0 or hi<0)})
    ctx = pd.DataFrame(context_rows)
    ctx.to_csv(ROOT / "04_model_context/model_context_effects.tsv", sep="\t", index=False)

    # Rank-score sensitivity for the same effects.
    rank_rows = []
    for sid in rank_score.columns:
        a, b = rank_score.loc[ng_ids,sid].to_numpy(), rank_score.loc[tr_ids,sid].to_numpy()
        diff, lo, hi = boot_diff(a,b,rng)
        rank_rows.append({"signature_id": sid, "median_difference_NextGen_minus_traditional": diff,
                          "bootstrap_CI_low": lo, "bootstrap_CI_high": hi, "hedges_g": hedges_g(a,b),
                          "direction_matches_z_score": np.sign(diff)==np.sign(ctx.loc[ctx.analysis.eq(f'score_z__{sid}'),'median_difference_NextGen_minus_traditional'].iloc[0])})
    pd.DataFrame(rank_rows).to_csv(ROOT / "11_sensitivity/rank_score_context_effects.tsv", sep="\t", index=False)

    # Gene-by-context table, patient frozen reference and candidate culture enrichment under frozen thresholds.
    pr = pref.set_index("gene")
    patient_primary = pr["GBmap_malignant_pseudobulk"].reindex(mes_avail).astype(float)
    patient_sens = pr["GSE174554_tumor_pseudobulk"].reindex(mes_avail).astype(float)
    patient_z = (patient_primary-patient_primary.mean())/patient_primary.std(ddof=1)
    patient_sens_z = (patient_sens-patient_sens.mean())/patient_sens.std(ddof=1)
    g_rows = []
    for g in mes_avail:
        a, b = z.loc[ng_ids,mapped_col[g]].to_numpy(float), z.loc[tr_ids,mapped_col[g]].to_numpy(float)
        diff, lo, hi = boot_diff(a,b,rng,fun=np.mean)
        rank_p = rankdata(np.r_[patient_z.values, [np.nanmean(a), np.nanmean(b)]], nan_policy="omit")
        # Cross-context culture candidate: traditional > NextGen, stable CI, large effect, and traditional rank exceeds both references.
        vals3 = pd.Series({"patient": patient_z.get(g,np.nan), "NextGen": np.nanmean(a), "traditional": np.nanmean(b)})
        pct = vals3.rank(pct=True)
        culture = bool(hedges_g(b,a) >= .8 and lo < 0 and hi < 0 and pct["traditional"]-max(pct["patient"],pct["NextGen"]) >= .25)
        g_rows.append({"gene": g, "patient_reference_primary": patient_primary.get(g,np.nan), "patient_reference_z": patient_z.get(g,np.nan),
                       "patient_reference_sensitivity": patient_sens.get(g,np.nan), "NextGen_mean_z": np.nanmean(a), "traditional_mean_z": np.nanmean(b),
                       "NextGen_median_z": np.nanmedian(a), "traditional_median_z": np.nanmedian(b),
                       "mean_difference_NextGen_minus_traditional": diff, "bootstrap_CI_low": lo, "bootstrap_CI_high": hi,
                       "hedges_g_NextGen_minus_traditional": hedges_g(a,b), "provenance_class": prov.get(g,"UNRESOLVED"),
                       "modality_dependent": g in mod_genes, "culture_enriched_candidate": culture,
                       "culture_candidate_rule": "g(traditional-NextGen)>=0.8; CI traditional>NextGen; traditional relative rank exceeds patient and NextGen by >=0.25"})
    gene_ctx = pd.DataFrame(g_rows)
    gene_ctx.to_csv(ROOT / "05_mes95_contribution/mes95_gene_by_context.tsv", sep="\t", index=False)

    cont = mes_z.loc[mids].copy()
    cont.insert(0,"NextGen_or_traditional",reg.loc[mids,"NextGen_or_traditional"])
    cont.reset_index(names="ModelID").to_csv(ROOT / "05_mes95_contribution/mes95_contribution_matrix.tsv", sep="\t", index=False)
    gene_ctx.to_csv(ROOT / "05_mes95_contribution/mes95_contribution_shift.tsv", sep="\t", index=False)
    comp.loc[mids].groupby(reg.loc[mids,"NextGen_or_traditional"]).agg(["mean","median"]).to_csv(ROOT / "05_mes95_contribution/provenance_by_context.tsv", sep="\t")

    # Patient concordance and direct test of patient-NextGen versus patient-traditional correlation.
    def group_corr(ids0, patient_vec):
        v = mes_z.loc[ids0].mean(axis=0).reindex(mes_avail).to_numpy()
        return corr_safe(patient_vec.to_numpy(),v,"spearman"), corr_safe(patient_vec.to_numpy(),v,"pearson")
    ng_sp, ng_pe = group_corr(ng_ids,patient_z); tr_sp,tr_pe = group_corr(tr_ids,patient_z)
    boots = np.empty(NBOOT)
    for i in range(NBOOT):
        nmean = mes_z.loc[rng.choice(ng_ids,len(ng_ids),replace=True)].mean(axis=0).to_numpy()
        tmean = mes_z.loc[rng.choice(tr_ids,len(tr_ids),replace=True)].mean(axis=0).to_numpy()
        boots[i] = corr_safe(patient_z.to_numpy(),nmean)-corr_safe(patient_z.to_numpy(),tmean)
    pooled_ids = np.array(ng_ids+tr_ids)
    perm = np.empty(NPERM)
    for i in range(NPERM):
        sh = rng.permutation(pooled_ids)
        perm[i] = corr_safe(patient_z.to_numpy(),mes_z.loc[sh[:len(ng_ids)]].mean().to_numpy())-corr_safe(patient_z.to_numpy(),mes_z.loc[sh[len(ng_ids):]].mean().to_numpy())
    conc = pd.DataFrame([
        {"patient_reference":"GBmap_malignant_pseudobulk","context":"NextGen","spearman":ng_sp,"pearson":ng_pe,"n_models":len(ng_ids)},
        {"patient_reference":"GBmap_malignant_pseudobulk","context":"traditional","spearman":tr_sp,"pearson":tr_pe,"n_models":len(tr_ids)},
    ])
    conc.to_csv(ROOT / "05_mes95_contribution/patient_context_concordance.tsv",sep="\t",index=False)
    pd.DataFrame([{"contrast":"patient_NextGen_minus_patient_traditional_Spearman","observed_delta":ng_sp-tr_sp,
                   "bootstrap_CI_low":np.quantile(boots,.025),"bootstrap_CI_high":np.quantile(boots,.975),
                   "permutation_p_two_sided":(1+np.sum(np.abs(perm)>=abs(ng_sp-tr_sp)))/(NPERM+1),
                   "patient_reference_scope":"frozen malignant-pseudobulk gene vector; not complete bulk tumour"}])\
        .to_csv(ROOT / "05_mes95_contribution/patient_nextgen_vs_traditional_concordance.tsv",sep="\t",index=False)

    # Two prespecified same-score matching strategies.
    scores95 = score["CANONICAL_MES95"]
    A, B = scores95.loc[ng_ids].to_numpy(), scores95.loc[tr_ids].to_numpy()
    ri, ci = linear_sum_assignment(np.abs(A[:,None]-B[None,:]))
    nearest_pairs = [(ng_ids[i],tr_ids[j]) for i,j in zip(ri,ci)]
    cal = .2*np.std(np.r_[A,B],ddof=1)
    cal_idx = greedy_caliper(A,B,cal)
    cal_pairs = [(ng_ids[i],tr_ids[j]) for i,j in cal_idx]
    contribution = mes_z.loc[mids,mes_avail]
    components = comp.loc[mids,exclusive]
    pair_df = pd.concat([
        pair_metrics([x[0] for x in nearest_pairs],[x[1] for x in nearest_pairs],contribution,components,scores95,"nearest_optimal"),
        pair_metrics([x[0] for x in cal_pairs],[x[1] for x in cal_pairs],contribution,components,scores95,"caliper_0.2SD")
    ],ignore_index=True)
    pair_df.to_csv(ROOT / "06_matched_mes_composition/matched_same_mes_composition.tsv",sep="\t",index=False)

    # Label-permutation null reruns score matching and preserves joint score/contribution structure.
    null_rows=[]
    all_ids=np.array(ng_ids+tr_ids)
    nng=len(ng_ids)
    for it in range(1000):
        sh=rng.permutation(all_ids); aa=sh[:nng].tolist(); bb=sh[nng:].tolist()
        av=scores95.loc[aa].to_numpy(); bv=scores95.loc[bb].to_numpy()
        ii,jj=linear_sum_assignment(np.abs(av[:,None]-bv[None,:]))
        tmp=pair_metrics([aa[i] for i in ii],[bb[j] for j in jj],contribution,components,scores95,"null")
        null_rows.append({"iteration":it+1,"mean_contribution_spearman":tmp.contribution_spearman.mean(),
                          "mean_cosine_similarity":tmp.contribution_cosine_similarity.mean(),
                          "mean_rank_distance":tmp.contribution_rank_distance.mean(),
                          "mean_provenance_distance":tmp.provenance_euclidean_distance.mean()})
    null=pd.DataFrame(null_rows); null.to_csv(ROOT / "06_matched_mes_composition/matched_null.tsv",sep="\t",index=False)
    obs=pair_df.loc[pair_df.method.eq("nearest_optimal")]
    summ=[]
    for oc,nc,direction in [("contribution_spearman","mean_contribution_spearman","lower"),("contribution_cosine_similarity","mean_cosine_similarity","lower"),
                            ("contribution_rank_distance","mean_rank_distance","higher"),("provenance_euclidean_distance","mean_provenance_distance","higher")]:
        v=obs[oc].mean(); nv=null[nc].to_numpy()
        p=(1+np.sum(nv<=v))/(len(nv)+1) if direction=="lower" else (1+np.sum(nv>=v))/(len(nv)+1)
        summ.append({"metric":oc,"observed_mean":v,"null_mean":nv.mean(),"null_CI_low":np.quantile(nv,.025),"null_CI_high":np.quantile(nv,.975),
                     "one_sided_permutation_p":p,"alternative":direction,"matched_pair_n":len(obs),"caliper":cal})
    pd.DataFrame(summ).to_csv(ROOT / "06_matched_mes_composition/matched_composition_summary.tsv",sep="\t",index=False)

    receipt={"seed":SEED,"bootstrap_n":NBOOT,"permutation_n":NPERM,"expression_scale":"exported log2(TPM+1); no repeat transform",
             "primary_scoring":"pooled gene-wise z-score then frozen weighted mean","strict_RNA":{"NextGen":len(ng_ids),"traditional":len(tr_ids)},
             "same_score_caliper":cal,"patient_reference":"GBmap_malignant_pseudobulk"}
    (ROOT/"17_logs/phase2_expression_receipt.json").write_text(json.dumps(receipt,indent=2),encoding="utf-8")


if __name__ == "__main__":
    main()
