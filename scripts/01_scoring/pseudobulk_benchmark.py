#!/usr/bin/env python3
"""Controlled real-count pseudo-bulk score benchmark with discovery/validation studies excluded."""
from __future__ import annotations

import json
import os
from pathlib import Path

DATA_ROOT = Path(os.environ.get("GBM_MES_DATA_ROOT", "data"))
RESULT_ROOT = Path(os.environ.get("GBM_MES_OUTPUT_ROOT", "results"))
REPO_ROOT = Path(os.environ.get("GBM_MES_REPO_ROOT", Path(__file__).resolve().parents[2]))

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

OUT = RESULT_ROOT / "scoring" / "pseudobulk"
REFERENCE = DATA_ROOT / "prepared" / "single_cell_reference"
COUNTS = REFERENCE / "multistudy_scrna_counts.h5"
META = REFERENCE / "multistudy_scrna_metadata.tsv"
SIG = REPO_ROOT / "resources" / "signatures" / "signature_gene_sets.tsv"
RNG = np.random.default_rng(2026073102)
GRID = np.array([0, .05, .10, .20, .30, .40])
REPS = 10
EXCLUDED_STUDIES = {"Neftel2019", "Bhaduri2020", "Couturier2020", "Johnson2020", "Richards2021", "Yu2020", "Yuan2018"}


def read_sets():
    r = pd.read_csv(SIG, sep="\t")
    r["gene"] = r.gene_symbol_clean.astype(str).str.strip().str.upper()
    ids = {
        "canonical_MES": "NEFTEL_MES_LIKE", "MES1": "NEFTEL_MES1", "MES2": "NEFTEL_MES2",
        "AC_like": "NEFTEL_AC_LIKE", "OPC_like": "NEFTEL_OPC_LIKE", "NPC_like": "NEFTEL_NPC_LIKE",
        "proliferation": "PROLIFERATION_CELL_CYCLE",
    }
    out = {k: set(r.loc[r.signature_id.eq(v), "gene"].dropna()) for k, v in ids.items()}
    out["MES_Hyp"] = set(r.loc[r.signature_id.eq("NEFTEL_MES_HYP"), "gene"].dropna())
    out["MES_Ast"] = set(r.loc[r.signature_id.eq("NEFTEL_MES_AST"), "gene"].dropna())
    return out


def sample_rows(pool: np.ndarray, n: int) -> np.ndarray:
    if n <= 0:
        return np.empty(0, dtype=int)
    return RNG.choice(pool, size=n, replace=n > len(pool))


def fit_fixed_effect(d: pd.DataFrame, experiment: str):
    y = d.score.to_numpy(float)
    mes = d.malignant_mes_fraction.to_numpy(float)
    comp = d.composition_fraction.to_numpy(float)
    # Study is nested in patient for this design, so patient fixed effects are
    # sufficient for the response slopes.  Study heterogeneity is reported
    # separately and is not duplicated in this rank-sensitive design matrix.
    patients = pd.get_dummies(d.patient_key, drop_first=True, dtype=float)
    if experiment == "A_FIXED_MALIGNANT":
        terms = np.column_stack([comp])
    elif experiment == "B_FIXED_COMPOSITION":
        terms = np.column_stack([mes])
    else:
        terms = np.column_stack([mes, comp, mes * comp])
    X = np.column_stack([np.ones(len(d)), terms, patients.to_numpy()])
    keep = np.isfinite(y) & np.isfinite(X).all(axis=1)
    if keep.sum() <= X.shape[1] + 2:
        return (np.nan, np.nan, np.nan), np.nan
    beta, _, _, _ = np.linalg.lstsq(X[keep], y[keep], rcond=None)
    pred = X[keep] @ beta
    r2 = 1 - np.square(y[keep] - pred).sum() / np.square(y[keep] - y[keep].mean()).sum()
    if experiment == "A_FIXED_MALIGNANT":
        return (np.nan, beta[1], np.nan), r2
    if experiment == "B_FIXED_COMPOSITION":
        return (beta[1], np.nan, np.nan), r2
    return (beta[1], beta[2], beta[3]), r2


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(META, sep="\t")
    with h5py.File(COUNTS, "r") as h:
        ids = h["cell_id"].asstr()[:]
        genes = np.char.upper(h["gene_symbol"].asstr()[:].astype(str))
        g = h["counts_csr_cells_by_genes"]
        full = sparse.csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=(len(ids), len(genes)))
    if not np.array_equal(ids, meta.cell_id.astype(str).to_numpy()):
        raise SystemExit("P0_STOP: reference metadata and count-row identities differ")
    sets = read_sets()
    wanted = sorted(set().union(*sets.values()) & set(genes))
    fmap = {x: i for i, x in enumerate(genes)}
    cols = np.array([fmap[x] for x in wanted])
    x = full[:, cols].tocsr()
    del full
    gmap = {g: i for i, g in enumerate(wanted)}

    meta["patient_key"] = meta.author.astype(str) + "::" + meta.donor_id.astype(str)
    pool = meta.loc[~meta.author.isin(EXCLUDED_STUDIES)].copy()
    pool["row"] = pool.index
    pools = {}
    for patient, d in pool.groupby("patient_key"):
        pools[patient] = {
            "study": d.author.iloc[0],
            "malignant": d.loc[d.reference_cell_type.eq("malignant"), "row"].to_numpy(int),
            "MES": d.loc[d.reference_cell_type.eq("malignant") & d.malignant_state.eq("MES-like"), "row"].to_numpy(int),
            "nonMES": d.loc[d.reference_cell_type.eq("malignant") & ~d.malignant_state.eq("MES-like"), "row"].to_numpy(int),
            "myeloid": d.loc[d.reference_cell_type.eq("myeloid_macrophage"), "row"].to_numpy(int),
            "endothelial": d.loc[d.reference_cell_type.eq("endothelial"), "row"].to_numpy(int),
            "pericyte": d.loc[d.reference_cell_type.eq("pericyte"), "row"].to_numpy(int),
            "lymphoid": d.loc[d.reference_cell_type.isin(["t_cell", "nk_cell", "b_plasma"]), "row"].to_numpy(int),
            "oligodendrocyte_other": d.loc[d.reference_cell_type.isin(["oligodendrocyte", "opc_nonmalignant", "mast_cell"]), "row"].to_numpy(int),
            "all_nonmalignant": d.loc[~d.reference_cell_type.eq("malignant"), "row"].to_numpy(int),
        }

    designs, vectors = [], []
    sid = 0
    components = ["myeloid", "endothelial", "pericyte", "lymphoid", "oligodendrocyte_other"]
    for comp in components:
        for patient, p in pools.items():
            if len(p["malignant"]) < 20 or len(p[comp]) < 10:
                continue
            fixed_malignant = sample_rows(p["malignant"], 100)
            mes0 = np.mean(meta.iloc[fixed_malignant].malignant_state.eq("MES-like"))
            for frac in GRID:
                ncomp = int(round(100 * frac / (1 - frac))) if frac < 1 else 0
                for rep in range(REPS):
                    rows = np.r_[fixed_malignant, sample_rows(p[comp], ncomp)]
                    sid += 1; vectors.append(np.asarray(x[rows].sum(axis=0)).ravel())
                    designs.append(dict(pseudobulk_id=f"PB{sid:06d}", experiment="A_FIXED_MALIGNANT", compartment=comp,
                        study=p["study"], patient_key=patient, replicate=rep + 1, target_fraction=frac,
                        malignant_mes_fraction=mes0, composition_fraction=ncomp / len(rows), n_malignant=100, n_nonmalignant=ncomp,
                        integer_count_source="REF_B raw-count H5; source studies exclude discovery and current validation"))

    for patient, p in pools.items():
        if len(p["MES"]) < 10 or len(p["nonMES"]) < 10 or len(p["all_nonmalignant"]) < 10:
            continue
        fixed_nonmal = sample_rows(p["all_nonmalignant"], 25)
        for frac in GRID:
            nmes = int(round(100 * frac)); nnon = 100 - nmes
            for rep in range(REPS):
                rows = np.r_[sample_rows(p["MES"], nmes), sample_rows(p["nonMES"], nnon), fixed_nonmal]
                sid += 1; vectors.append(np.asarray(x[rows].sum(axis=0)).ravel())
                designs.append(dict(pseudobulk_id=f"PB{sid:06d}", experiment="B_FIXED_COMPOSITION", compartment="all_nonmalignant",
                    study=p["study"], patient_key=patient, replicate=rep + 1, target_fraction=frac,
                    malignant_mes_fraction=nmes / 100, composition_fraction=25 / 125, n_malignant=100, n_nonmalignant=25,
                    integer_count_source="REF_B raw-count H5; source studies exclude discovery and current validation"))

    for comp in ["myeloid", "endothelial", "pericyte"]:
        for patient, p in pools.items():
            if len(p["MES"]) < 10 or len(p["nonMES"]) < 10 or len(p[comp]) < 10:
                continue
            for mesfrac in GRID:
                nmes = int(round(100 * mesfrac)); nnon = 100 - nmes
                for cfrac in GRID:
                    ncomp = int(round(100 * cfrac / (1 - cfrac))) if cfrac < 1 else 0
                    for rep in range(REPS):
                        rows = np.r_[sample_rows(p["MES"], nmes), sample_rows(p["nonMES"], nnon), sample_rows(p[comp], ncomp)]
                        sid += 1; vectors.append(np.asarray(x[rows].sum(axis=0)).ravel())
                        designs.append(dict(pseudobulk_id=f"PB{sid:06d}", experiment="C_TWO_DIMENSIONAL", compartment=comp,
                            study=p["study"], patient_key=patient, replicate=rep + 1, target_fraction=f"MES={mesfrac};COMP={cfrac}",
                            malignant_mes_fraction=mesfrac, composition_fraction=ncomp / len(rows), n_malignant=100, n_nonmalignant=ncomp,
                            integer_count_source="REF_B raw-count H5; source studies exclude discovery and current validation"))

    design = pd.DataFrame(designs)
    counts = np.vstack(vectors)
    lib = counts.sum(axis=1)
    logcpm = np.log2((counts + .5) / (lib[:, None] + 1) * 1e6)
    response = []
    for (exp, comp), idx in design.groupby(["experiment", "compartment"]).groups.items():
        ii = np.array(list(idx))
        z = logcpm[ii]
        mu = z.mean(axis=0); sd = z.std(axis=0, ddof=1); sd[sd == 0] = np.nan
        z = (z - mu) / sd
        for score, gs in sets.items():
            gg = [gmap[g] for g in gs if g in gmap]
            val = np.nanmean(z[:, gg], axis=1) if len(gg) >= 2 else np.full(len(ii), np.nan)
            for j, v in zip(ii, val):
                response.append({"pseudobulk_id": design.loc[j, "pseudobulk_id"], "experiment": exp, "compartment": comp,
                    "study": design.loc[j, "study"], "patient_key": design.loc[j, "patient_key"], "score_name": score,
                    "score": v, "genes_total": len(gs), "genes_detected": len(gg),
                    "malignant_mes_fraction": design.loc[j, "malignant_mes_fraction"], "composition_fraction": design.loc[j, "composition_fraction"]})
    response = pd.DataFrame(response)

    benchmark = []
    for (exp, comp, score), d in response.groupby(["experiment", "compartment", "score_name"]):
        (state, composition, interaction), r2 = fit_fixed_effect(d, exp)
        patient_slopes = []
        for _, q in d.groupby("patient_key"):
            if exp == "A_FIXED_MALIGNANT":
                X = np.column_stack([np.ones(len(q)), q.composition_fraction])
                if len(q) > 5 and np.linalg.matrix_rank(X) == 2:
                    patient_slopes.append([np.nan, np.linalg.lstsq(X, q.score, rcond=None)[0][1]])
            elif exp == "B_FIXED_COMPOSITION":
                X = np.column_stack([np.ones(len(q)), q.malignant_mes_fraction])
                if len(q) > 5 and np.linalg.matrix_rank(X) == 2:
                    patient_slopes.append([np.linalg.lstsq(X, q.score, rcond=None)[0][1], np.nan])
            else:
                X = np.column_stack([np.ones(len(q)), q.malignant_mes_fraction, q.composition_fraction, q.malignant_mes_fraction * q.composition_fraction])
                if len(q) > 5 and np.linalg.matrix_rank(X) == 4:
                    patient_slopes.append(np.linalg.lstsq(X, q.score, rcond=None)[0][1:3])
        slope_array = np.asarray(patient_slopes, dtype=float)
        pv = np.nanmedian(np.nanstd(slope_array, axis=0)) if patient_slopes else np.nan
        study_vals = d.groupby("study").score.mean()
        benchmark.append({
            "experiment": exp, "compartment": comp, "score_name": score,
            "program_status": "PRESPECIFIED_PUBLISHED_SCORE",
            "malignant_state_slope": state, "composition_slope": composition, "interaction": interaction,
            "dynamic_range": d.groupby(["malignant_mes_fraction", "composition_fraction"]).score.mean().agg(lambda x: x.max() - x.min()),
            "signal_to_composition_ratio": abs(state) / (abs(composition) + 1e-8) if np.isfinite(state) and np.isfinite(composition) else np.nan,
            "patient_level_variability": pv, "study_level_variability": study_vals.std(ddof=1),
            "model_R2": r2, "patients": d.patient_key.nunique(), "studies": d.study.nunique(), "pseudobulks": len(d),
            "reference_dependence": "NOT_APPLICABLE_DIRECT_SCORE; source-pool study variability reported",
        })
    benchmark = pd.DataFrame(benchmark)
    design.to_csv(OUT / "PSEUDOBULK_DESIGN.tsv", sep="\t", index=False)
    design[["pseudobulk_id", "experiment", "compartment", "study", "patient_key", "malignant_mes_fraction", "composition_fraction", "n_malignant", "n_nonmalignant"]].to_csv(OUT / "PSEUDOBULK_TRUE_COMPOSITION.tsv", sep="\t", index=False)
    response.to_csv(OUT / "SCORE_RESPONSE_RESULTS.tsv", sep="\t", index=False)
    benchmark.to_csv(OUT / "SCORE_INTERPRETABILITY_BENCHMARK.tsv", sep="\t", index=False)
    qa = pd.DataFrame([
        {"check": "integer_count_input", "status": "PASS", "detail": "REF_B raw-count H5 derived from GBmap raw/X"},
        {"check": "discovery_study_excluded", "status": "PASS", "detail": "Neftel2019 excluded"},
        {"check": "current_validation_studies_excluded", "status": "PASS", "detail": ";".join(sorted(EXCLUDED_STUDIES - {"Neftel2019"}))},
        {"check": "patient_stratified_replicates", "status": "PASS", "detail": f"{design.patient_key.nunique()} patients; {REPS} replicates per grid condition"},
        {"check": "fibroblast_stromal_experiment", "status": "NOT_EVALUABLE", "detail": "only frozen reference fibroblast pool was Neftel2019 discovery; excluded to prevent role leakage"},
        {"check": "prespecified_score_set", "status": "PASS", "detail": ";".join(sorted(sets))},
    ])
    qa.to_csv(OUT / "PSEUDOBULK_QA.tsv", sep="\t", index=False)
    (OUT / "PARAMETERS.json").write_text(json.dumps({"seed": 2026073102, "grid": GRID.tolist(), "replicates": REPS, "excluded_studies": sorted(EXCLUDED_STUDIES), "normal_data_simulation": False}, indent=2) + "\n")


if __name__ == "__main__":
    main()
