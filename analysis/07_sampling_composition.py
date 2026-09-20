"""Frozen CARE sampling-null and non-additive composition decomposition.
Source-dependent execution: external frozen project layout; no source rows shipped.
Outputs include patient-level intermediates: write to an external private directory.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
EPS=1e-12
RNG=np.random.default_rng(202608244)
BOOT=1000
SPECS=(SimpleNamespace(dataset_id="DS01",label="CARE",scope="ALL_AUTHOR_ANNOTATED",membership="AS_CARE_LOPO56.tsv",expected_n=56),)
def orient_loadings(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, float).copy()
    for j in range(v.shape[1]):
        idx = int(np.argmax(np.abs(v[:, j])))
        if v[idx, j] < 0:
            v[:, j] *= -1
    return v

def fit_reference(x: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    _, s, vt = np.linalg.svd(x, full_matrices=False)
    v = orient_loadings(vt[:k].T)
    all_eig = (s**2) / max(x.shape[0] - 1, 1)
    eig = all_eig[:k]
    explained = eig / all_eig.sum() if all_eig.sum() > 0 else np.full(k, np.nan)
    return v, eig, explained

def read_pair_matrix(spec: DatasetSpec) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    scores = pd.read_parquet(V1 / f"04_model_scores/{spec.dataset_id}/PRIMARY201_MODEL_SCORES.parquet")
    scores = scores[scores["score_scope"].eq(spec.scope)].copy()
    membership = pd.read_csv(V2 / "00_governance/freeze/memberships" / spec.membership, sep="\t")
    if "primary_sample_id" not in membership:
        pair = pd.read_csv(V1 / "05_modules/module_C_longitudinal/PATIENT_LONGITUDINAL_SIGNATURE_SIMILARITY.tsv", sep="\t")
        pair = pair[(pair["dataset_id"].eq(spec.dataset_id)) & pair["score_scope"].eq(spec.scope)][
            ["patient_id", "primary_sample_id", "recurrent_sample_id"]
        ].drop_duplicates()
        membership = membership.merge(pair, on="patient_id", how="left", validate="one_to_one")
    models = pd.read_csv(V2 / "00_governance/freeze/memberships/AS_PRIMARY201_V2.tsv", sep="\t")[
        "model_instance_id"
    ].astype(str).tolist()
    index = scores.set_index(["patient_id", "sample_id", "model_instance_id"])["model_score"]
    p_rows, r_rows, patients = [], [], []
    for row in membership.itertuples(index=False):
        pid = str(row.patient_id)
        ps, rs = str(row.primary_sample_id), str(row.recurrent_sample_id)
        try:
            pv = index.loc[(pid, ps)].reindex(models)
            rv = index.loc[(pid, rs)].reindex(models)
        except KeyError:
            continue
        p_rows.append(pv.to_numpy(float))
        r_rows.append(rv.to_numpy(float))
        patients.append(pid)
    xp = pd.DataFrame(p_rows, index=patients, columns=models)
    xr = pd.DataFrame(r_rows, index=patients, columns=models)
    complete = xp.notna().all(axis=0) & xr.notna().all(axis=0)
    xp, xr = xp.loc[:, complete], xr.loc[:, complete]
    if len(xp) != spec.expected_n or not xp.index.equals(xr.index) or xp.shape[1] < 200:
        raise RuntimeError(f"{spec.label} pair/matrix identity failure: {xp.shape}, {xr.shape}")
    return xp, xr, xp.columns.tolist()
def reference_for_patient(xp: np.ndarray, i: int, k: int):
    train = xp[np.arange(len(xp)) != i]
    mu = train.mean(0); sd = train.std(0, ddof=1)
    valid = np.isfinite(mu) & np.isfinite(sd) & (sd > EPS)
    z = (train[:, valid] - mu[valid]) / sd[valid]
    v, _, _ = fit_reference(z, k)
    return mu, sd, valid, v

def decompose(delta: np.ndarray, v: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    parallel = delta @ v @ v.T
    orth = delta - parallel
    total = np.sum(delta * delta, axis=-1)
    pe = np.sum(parallel * parallel, axis=-1)
    oe = np.sum(orth * orth, axis=-1)
    return pe, oe, total

def bootstrap_median(x: np.ndarray):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    draws = np.array([np.median(x[RNG.integers(0, len(x), len(x))]) for _ in range(BOOT)])
    return float(np.median(x)), float(np.quantile(draws, .025)), float(np.quantile(draws, .975))

def build_model_weights(gene_order: list[str]):
    models = pd.read_csv(V1 / "00_governance/frozen_parent_inputs/01_PRIMARY201_MODEL_MANIFEST_V2_1_1.tsv", sep="\t")
    mg = pd.read_csv(V1 / "00_governance/frozen_parent_inputs/02_PRIMARY201_GENE_LONG_V2_1_1.tsv", sep="\t")
    cw = pd.read_csv(V1 / "00_governance/frozen_parent_inputs/PRIMARY201_SOURCE_TO_CANONICAL_IDENTITY_CROSSWALK_V2_1_3.tsv", sep="\t")
    mg["source_scoring_feature_id"] = mg.scoring_feature_id.astype(str)
    cw["source_scoring_feature_id"] = cw.source_scoring_feature_id.astype(str)
    m = mg.merge(cw[["source_scoring_feature_id", "canonical_gene_symbol"]], on="source_scoring_feature_id", how="left")
    m["gene"] = m.canonical_gene_symbol.fillna(m.source_scoring_feature_id).str.upper()
    order = models.model_instance_id.astype(str).tolist(); gi = {g:i for i,g in enumerate(gene_order)}
    w = np.zeros((len(gene_order), len(order)))
    for j, mid in enumerate(order):
        genes = sorted(set(m.loc[m.model_instance_id == mid, "gene"]))
        present = [g for g in genes if g in gi]
        if len(present) >= 3 and len(present)/len(genes) >= .60:
            w[[gi[g] for g in present], j] = 1/len(present)
    return order, w
def main():
    sample_dir = OUTPUT / "03_sampling_null"; comp_dir = OUTPUT / "05_composition_intrinsic"
    sample_dir.mkdir(parents=True,exist_ok=True); comp_dir.mkdir(parents=True,exist_ok=True)
    xp_df, xr_df, model_order = read_pair_matrix(SPECS[0])
    xp, xr = xp_df.to_numpy(float), xr_df.to_numpy(float)
    pairs = pd.read_csv(V1 / "05_modules/module_C_longitudinal/PATIENT_LONGITUDINAL_SIGNATURE_SIMILARITY.tsv", sep="\t")
    pairs = pairs[(pairs.dataset_id == "DS01") & (pairs.score_scope == "ALL_AUTHOR_ANNOTATED")][
        ["patient_id", "primary_sample_id", "recurrent_sample_id"]].drop_duplicates().set_index("patient_id").loc[xp_df.index]

    split = pd.read_parquet(V3 / "03_sampling_null" / "CARE_SPLIT_HALF_MODEL_SCORES.parquet")
    down = pd.read_parquet(V3 / "03_sampling_null" / "CARE_MATCHED_CELL_DOWNSAMPLE_MODEL_SCORES.parquet")
    null_rows, compare_rows, down_rows = [], [], []
    for i, pid in enumerate(xp_df.index):
        pair = pairs.loc[pid]
        for k in [2, 7]:
            mu, sd, valid, v = reference_for_patient(xp, i, k)
            real_delta = (xr[i, valid] - xp[i, valid]) / sd[valid]
            rpe, roe, rte = decompose(real_delta[None, :], v)
            pooled = []
            for phase, sid in [("PRIMARY", pair.primary_sample_id), ("RECURRENT", pair.recurrent_sample_id)]:
                ss = split[split.sample_id == sid].sort_values(["replicate", "half"])
                a = ss[ss.half == "A"].sort_values("replicate")[np.asarray(model_order)[valid]].to_numpy(float)
                b = ss[ss.half == "B"].sort_values("replicate")[np.asarray(model_order)[valid]].to_numpy(float)
                d = (b-a)/sd[valid]
                pe, oe, te = decompose(d, v)
                pooled.append(np.sqrt(oe))
                for rep, p0, o0, t0 in zip(range(1, len(pe)+1), pe, oe, te):
                    null_rows.append({"patient_id":pid,"reference_k":k,"phase":phase,"replicate":rep,
                                      "parallel_energy":p0,"orthogonal_energy":o0,"total_energy":t0,
                                      "residual_magnitude":np.sqrt(o0),
                                      "orthogonal_energy_fraction":o0/t0 if t0>EPS else np.nan})
            null = np.concatenate(pooled)
            compare_rows.append({
                "patient_id":pid,"reference_k":k,"real_residual_magnitude":np.sqrt(roe[0]),
                "real_offspace_fraction":roe[0]/rte[0],"sampling_null_median_residual_magnitude":np.median(null),
                "sampling_null_95_high_residual_magnitude":np.quantile(null,.95),
                "biological_deformation_excess":np.sqrt(roe[0])-np.median(null),
                "standardized_excess":(np.sqrt(roe[0])-np.mean(null))/np.std(null,ddof=1),
                "empirical_p_real_le_sampling_null":(1+np.sum(null>=np.sqrt(roe[0])))/(1+len(null)),
                "null_replicate_n":len(null),"inference_unit":"patient"})

            pdn = down[down.sample_id == pair.primary_sample_id].sort_values("replicate")
            rdn = down[down.sample_id == pair.recurrent_sample_id].sort_values("replicate")
            dr = (rdn[np.asarray(model_order)[valid]].to_numpy(float)-pdn[np.asarray(model_order)[valid]].to_numpy(float))/sd[valid]
            pe, oe, te = decompose(dr, v)
            for rep,p0,o0,t0 in zip(range(1,len(pe)+1),pe,oe,te):
                down_rows.append({"patient_id":pid,"reference_k":k,"replicate":rep,"parallel_energy":p0,
                                  "orthogonal_energy":o0,"total_energy":t0,"residual_magnitude":np.sqrt(o0),
                                  "orthogonal_energy_fraction":o0/t0 if t0>EPS else np.nan,
                                  "original_residual_magnitude":np.sqrt(roe[0]),
                                  "original_offspace_fraction":roe[0]/rte[0]})
    null_df=pd.DataFrame(null_rows); compare=pd.DataFrame(compare_rows); down_df=pd.DataFrame(down_rows)
    null_df.to_csv(sample_dir/"SPLIT_HALF_DEFORMATION_NULL.tsv",sep="\t",index=False)
    compare.to_csv(sample_dir/"REAL_VS_SAMPLING_DEFORMATION.tsv",sep="\t",index=False)
    down_df.to_csv(sample_dir/"MATCHED_CELL_DOWNSAMPLE_SENSITIVITY.tsv",sep="\t",index=False)
    cohort=[]
    for k,g in compare.groupby("reference_k"):
        for col in ["real_residual_magnitude","sampling_null_median_residual_magnitude","biological_deformation_excess","standardized_excess"]:
            med,lo,hi=bootstrap_median(g[col].to_numpy(float)); cohort.append({"reference_k":k,"endpoint":col,"median":med,"ci_low":lo,"ci_high":hi,"patient_n":len(g)})
    pd.DataFrame(cohort).to_csv(sample_dir/"WP2_COHORT_SUMMARY.tsv",sep="\t",index=False)

    # Broad-compartment exact two-factor Shapley reconstruction.
    profiles = pd.read_parquet(V3 / "05_composition_intrinsic" / "CARE_BROAD_COMPARTMENT_REQUIRED_GENE_COUNTS.parquet")
    gene_cols = [c for c in profiles.columns if c not in {"sample_id","broad_compartment","cell_n","library_size"}]
    pb = pd.read_parquet(V1 / "03_processed_data/DS01_CARE/pseudobulk/DS01_REQUIRED_GENE_PSEUDOBULK_COUNTS.parquet")
    pb = pb[pb.score_scope == "ALL_AUTHOR_ANNOTATED"].copy()
    pb["expression"] = np.log2((pb["count"].astype(float)+.5)/(pb["library_size"].astype(float)+1)*1e6)
    wide = pb.pivot(index="sample_id",columns="gene_symbol",values="expression").sort_index(axis=1)
    gene_cols = wide.columns.tolist()
    gmu=wide.mean().to_numpy(float); gsd=wide.std(ddof=1).to_numpy(float); gvalid=gsd>EPS
    ordered_models,w=build_model_weights(gene_cols)
    if ordered_models != model_order: raise RuntimeError("Primary201 model-order drift")

    def score_mix(pvec, state_counts, state_lib):
        # Evaluate every counterfactual at the same large synthetic cell count.
        # This preserves the mixture's CPM while preventing the +0.5 count
        # pseudocount from becoming artificially dominant for per-cell means.
        synthetic_cell_n = 1_000_000.0
        counts = (pvec @ state_counts) * synthetic_cell_n
        lib = float(pvec @ state_lib) * synthetic_cell_n
        expr=np.log2((counts+.5)/(lib+1)*1e6); z=(expr-gmu)/gsd; z[~gvalid]=0
        return z@w

    class_df=pd.read_csv(V1/"00_governance/frozen_parent_results/PRIMARY201_FULL_GENE_CLASSIFICATION.tsv",sep="\t")
    pub=pd.read_csv(V2/"00_governance/freeze/memberships/AS_PRIMARY201_V2.tsv",sep="\t")
    meta=class_df[["model_instance_id","final_class"]].merge(pub[["model_instance_id","publication_id"]],on="model_instance_id")
    shapley_rows=[]; model_rows=[]
    for i,pid in enumerate(xp_df.index):
        pair=pairs.loc[pid]
        gp=profiles[profiles.sample_id==pair.primary_sample_id].set_index("broad_compartment")
        gr=profiles[profiles.sample_id==pair.recurrent_sample_id].set_index("broad_compartment")
        for threshold in [50,25]:
            all_comps=sorted(set(gp.index).union(gr.index))
            comps=[c for c in all_comps if c in gp.index and c in gr.index and
                   gp.loc[c,"cell_n"]>=threshold and gr.loc[c,"cell_n"]>=threshold]
            if len(comps)<2: continue
            low=[c for c in all_comps if c not in comps]
            def state_arrays(g):
                counts=[]; cells=[]; libs=[]
                for c in comps:
                    counts.append(g.loc[c,gene_cols].to_numpy(float)); cells.append(float(g.loc[c,"cell_n"])); libs.append(float(g.loc[c,"library_size"]))
                present_low=[c for c in low if c in g.index]
                if present_low:
                    counts.append(g.loc[present_low,gene_cols].sum(axis=0).to_numpy(float))
                    cells.append(float(g.loc[present_low,"cell_n"].sum())); libs.append(float(g.loc[present_low,"library_size"].sum()))
                return np.vstack(counts),np.asarray(cells),np.asarray(libs)
            cntP,cellP,libP=state_arrays(gp); cntR,cellR,libR=state_arrays(gr)
            # Both states have the same explicit labels plus one pooled-low label.
            # If one side has no low cells, append a zero-mass placeholder.
            if len(cellP)!=len(cellR):
                if len(cellP)<len(cellR): cntP=np.vstack([cntP,np.zeros(len(gene_cols))]); cellP=np.append(cellP,0.); libP=np.append(libP,0.)
                else: cntR=np.vstack([cntR,np.zeros(len(gene_cols))]); cellR=np.append(cellR,0.); libR=np.append(libR,0.)
            pP=cellP/cellP.sum(); pR=cellR/cellR.sum()
            EP=np.divide(cntP,cellP[:,None],out=np.zeros_like(cntP),where=cellP[:,None]>0)
            ER=np.divide(cntR,cellR[:,None],out=np.zeros_like(cntR),where=cellR[:,None]>0)
            lP=np.divide(libP,cellP,out=np.zeros_like(libP),where=cellP>0)
            lR=np.divide(libR,cellR,out=np.zeros_like(libR),where=cellR>0)
            PP=score_mix(pP,EP,lP); RP=score_mix(pR,EP,lP); PR=score_mix(pP,ER,lR); RR=score_mix(pR,ER,lR)
            C=.5*((RP-PP)+(RR-PR)); I=.5*((PR-PP)+(RR-RP)); delta=RR-PP
            identity=float(np.max(np.abs(delta-C-I)))
            for k in [2,7]:
                mu,sd,valid,v=reference_for_patient(xp,i,k)
                cz=C[valid]/sd[valid]; iz=I[valid]/sd[valid]; dz=delta[valid]/sd[valid]
                ape,aoe,ate=decompose(dz[None,:],v); cpe,coe,cte=decompose(cz[None,:],v); ipe,ioe,ite=decompose(iz[None,:],v)
                actual=(xr[i,valid]-xp[i,valid])/sd[valid]
                cos_actual=float(dz@actual/(np.linalg.norm(dz)*np.linalg.norm(actual))) if np.linalg.norm(dz)*np.linalg.norm(actual)>EPS else np.nan
                cos_ci=float(cz@iz/(np.linalg.norm(cz)*np.linalg.norm(iz))) if np.linalg.norm(cz)*np.linalg.norm(iz)>EPS else np.nan
                shapley_rows.append({"patient_id":pid,"cell_threshold":threshold,"included_compartment_n":len(comps),
                    "included_compartments":";".join(comps),"pooled_below_threshold_compartments":";".join(low),
                    "reference_k":k,"shapley_identity_max_abs_error":identity,
                    "reconstruction_vs_original_delta_cosine":cos_actual,"composition_energy":cte[0],"intrinsic_energy":ite[0],
                    "composition_intrinsic_cross_term":2*float(cz@iz),"total_reconstructed_energy":ate[0],
                    "composition_energy_fraction_nonadditive":cte[0]/(cte[0]+ite[0]) if cte[0]+ite[0]>EPS else np.nan,
                    "composition_intrinsic_cosine":cos_ci,"composition_parallel_energy":cpe[0],"composition_offspace_energy":coe[0],
                    "intrinsic_parallel_energy":ipe[0],"intrinsic_offspace_energy":ioe[0],"total_parallel_energy":ape[0],"total_offspace_energy":aoe[0]})
            for j,mid in enumerate(model_order):
                den=abs(C[j])+abs(I[j])
                model_rows.append({"patient_id":pid,"cell_threshold":threshold,"model_instance_id":mid,
                    "composition_component":C[j],"intrinsic_component":I[j],"reconstructed_delta":delta[j],
                    "composition_absolute_share":abs(C[j])/den if den>EPS else np.nan})
    shap=pd.DataFrame(shapley_rows); modeld=pd.DataFrame(model_rows).merge(meta,on="model_instance_id",how="left",validate="many_to_one")
    shap.to_csv(comp_dir/"COMPOSITION_INTRINSIC_SHAPLEY.tsv",sep="\t",index=False)
    modeld.to_csv(comp_dir/"COMPARTMENT_CLASS_DECOMPOSITION.tsv",sep="\t",index=False)


    g=compare[compare.reference_k==7]
    q=shap[(shap.reference_k==7)&(shap.cell_threshold==50)]
    t=shap[(shap.reference_k==7)&(shap.cell_threshold==25)]
    summary=dict(patient_n=len(g),observed_median=float(g.real_residual_magnitude.median()),null_median=float(g.sampling_null_median_residual_magnitude.median()),above_own_95=int((g.real_residual_magnitude>g.sampling_null_95_high_residual_magnitude).sum()),composition_n=len(q),composition_energy=float(q.composition_energy_fraction_nonadditive.median()),reconstruction_cosine=float(q.reconstruction_vs_original_delta_cosine.median()),sensitivity_n=len(t),sensitivity_energy=float(t.composition_energy_fraction_nonadditive.median()))
    (OUTPUT/"sampling_composition_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary))

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("input_root",help="Externally held frozen project layout");p.add_argument("--output",required=True)
    a=p.parse_args();ROOT=Path(a.input_root).resolve();OUTPUT=Path(a.output).resolve()
    repo=Path(__file__).resolve().parents[1]
    if OUTPUT==repo or repo in OUTPUT.parents:raise ValueError("Patient-level outputs must be outside public repository")
    V1=ROOT/"GBM_SIGNATURE_CONTEXT_STRESS_TEST_EXTENSION_V1_FINAL"
    V2=ROOT/"GBM_SIGNATURE_CONTEXT_EXTENSION_DEEPENING_V2_FINAL"
    V3=ROOT/"GBM_SIGNATURE_CONTEXT_DEFORMATION_DEEPENING_V3"
    if OUTPUT in [V1,V2,V3] or any(x in OUTPUT.parents for x in [V1,V2,V3]):raise ValueError("Do not write into frozen authorities")
    main()
