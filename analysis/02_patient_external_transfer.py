"""Frozen ecology transfer from prepared, externally held arrays; no data redistributed."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import SplineTransformer
SEED=20260829
BOOT_N=1000

def grouped_linear_cv(s: np.ndarray, y: np.ndarray, groups: np.ndarray, spline: bool = False) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    predictions = np.full_like(y, np.nan, dtype=float)
    sse = np.full(len(s), np.nan)
    sst = np.full(len(s), np.nan)
    for patient in np.unique(groups):
        test = groups == patient
        train = ~test
        xtr, xte = s[train, None], s[test, None]
        if spline:
            transform = SplineTransformer(n_knots=4, degree=3, knots="quantile", include_bias=False)
            xtr = transform.fit_transform(xtr)
            xte = transform.transform(xte)
        fit = LinearRegression().fit(xtr, y[train])
        predictions[test] = fit.predict(xte)
        center = np.nanmean(y[train], axis=0)
        sse[test] = np.sum((y[test] - predictions[test]) ** 2, axis=1)
        sst[test] = np.sum((y[test] - center) ** 2, axis=1)
    r2 = 1 - np.nansum(sse) / np.nansum(sst)
    return float(r2), predictions, sse, sst

def cluster_bootstrap_global(contrib: pd.DataFrame, model_ids: list[str], patient_order: list[str]) -> tuple[float, float, float, float, np.ndarray]:
    sse = np.zeros((len(patient_order), len(model_ids)))
    sst = np.zeros_like(sse)
    for i, p in enumerate(patient_order):
        q = contrib[contrib.patient_id.eq(p)].set_index("model_instance_id")
        sse[i] = q.sse.reindex(model_ids).to_numpy(float)
        sst[i] = q.sst.reindex(model_ids).to_numpy(float)
    rng = np.random.default_rng(SEED + 101)
    med = np.full(BOOT_N, np.nan)
    for b in range(BOOT_N):
        idx = rng.choice(len(patient_order), size=len(patient_order), replace=True)
        r2 = 1 - np.sum(sse[idx], axis=0) / np.where(np.sum(sst[idx], axis=0) > 0, np.sum(sst[idx], axis=0), np.nan)
        med[b] = np.nanmedian(r2)
    return float(np.nanmedian(med)), float(np.nanquantile(med, .025)), float(np.nanquantile(med, .975)), float(np.mean(np.isfinite(med))), med

def ecology_cv(scores: pd.DataFrame, coverage: pd.DataFrame, meta: pd.DataFrame, y: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    groups = meta.patient_id.to_numpy(str)
    rows, contrib_rows = [], []
    evaluable = coverage.loc[coverage.evaluable, "model_instance_id"].tolist()
    for model in evaluable:
        s = scores[model].to_numpy(float)
        r2, pred, sse, sst = grouped_linear_cv(s, y, groups, spline=False)
        r2_spline, _, _, _ = grouped_linear_cv(s, y, groups, spline=True)
        axis_r2 = []
        for j in range(y.shape[1]):
            _, p1, e1, t1 = grouped_linear_cv(s, y[:, [j]], groups, spline=False)
            axis_r2.append(1 - np.nansum(e1) / np.nansum(t1))
        rows.append({
            "model_instance_id": model, "dataset_id": "DS03_GBMSPACE", "patient_n": len(np.unique(groups)),
            "sample_n": len(groups), "CV_multivariate_R2": r2,
            "RESIDUAL_BIOLOGY_FRACTION": 1 - r2,
            "CV_R2_MES": axis_r2[0], "CV_R2_hypoxia": axis_r2[1],
            "CV_R2_matrix": axis_r2[2], "CV_R2_myeloid": axis_r2[3],
            "CV_multivariate_R2_spline_sensitivity": r2_spline,
            "primary_model": "linear scalar score to four frozen axes",
            "CV_scheme": "leave-one-patient-out; patient-exclusive folds",
            "nonlinearity_sensitivity": "fixed 4-knot cubic spline; no endpoint-based selection",
        })
        for p in np.unique(groups):
            idx = groups == p
            contrib_rows.append({"patient_id": p, "model_instance_id": model, "sse": float(np.nansum(sse[idx])), "sst": float(np.nansum(sst[idx]))})
    table = pd.DataFrame(rows)
    contrib = pd.DataFrame(contrib_rows)
    med, lo, hi, valid, draws = cluster_bootstrap_global(contrib, evaluable, sorted(np.unique(groups)))
    global_summary = {
        "domain": "FOUR_FROZEN_ECOLOGY_AXES", "model_n": len(evaluable), "patient_n": len(np.unique(groups)),
        "median_CV_R2": float(table.CV_multivariate_R2.median()),
        "patient_bootstrap_median_CV_R2": med, "bootstrap_CI_low": lo, "bootstrap_CI_high": hi,
        "median_residual_fraction": float(table.RESIDUAL_BIOLOGY_FRACTION.median()),
        "residual_bootstrap_CI_low": float(1 - hi), "residual_bootstrap_CI_high": float(1 - lo),
        "bootstrap_n": BOOT_N, "bootstrap_valid_fraction": valid,
    }
    return table, contrib, global_summary

if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("input",help="NPZ: scores (n,m), endpoints (n,4), groups (n), model_ids (m)")
    p.add_argument("--output",required=True)
    a=p.parse_args(); z=np.load(a.input,allow_pickle=False)
    scores=pd.DataFrame(z["scores"],columns=z["model_ids"].astype(str))
    y=z["endpoints"]; groups=z["groups"].astype(str)
    if y.shape!=(len(scores),4) or len(groups)!=len(scores) or not np.isfinite(scores).all().all() or not np.isfinite(y).all():
        raise ValueError("Invalid prepared ecology arrays")
    coverage=pd.DataFrame({"model_instance_id":scores.columns,"evaluable":True})
    _,_,summary=ecology_cv(scores,coverage,pd.DataFrame({"patient_id":groups}),y)
    Path(a.output).write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary))
