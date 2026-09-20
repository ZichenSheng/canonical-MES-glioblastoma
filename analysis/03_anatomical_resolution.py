"""Frozen tumour-blocked association and leave-one-tumour-out resolution."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
Q1_NPERM=10000
BRANCHES={"ALL5":["LE","IT","CT","PAN","MVP"],"CORE3":["CT","PAN","MVP"],"PAN_MVP":["PAN","MVP"]}

def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    x = p[ok]
    if not len(x):
        return out
    order = np.argsort(x)
    ranked = x[order] * len(x) / np.arange(1, len(x) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    tmp = np.empty_like(x)
    tmp[order] = np.clip(ranked, 0, 1)
    out[ok] = tmp
    return out

def one_hot(values: np.ndarray, levels: list[str]) -> np.ndarray:
    return np.column_stack([(values == x).astype(float) for x in levels[1:]])

def orthonormal_basis(x: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    if x.size == 0:
        return np.empty((x.shape[0], 0))
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    rank = int(np.sum(s > tol * max(x.shape) * (s[0] if len(s) else 1.0)))
    return u[:, :rank]

def permute_within_tumor(labels: np.ndarray, tumors: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    out = labels.copy()
    for tumor in np.unique(tumors):
        idx = np.flatnonzero(tumors == tumor)
        out[idx] = rng.permutation(out[idx])
    return out

def q1_informativeness(y: np.ndarray, labels: np.ndarray, tumors: np.ndarray,
                       levels: list[str], rng: np.random.Generator) -> pd.DataFrame:
    tumor_levels = list(pd.unique(tumors))
    d = np.column_stack([np.ones(len(tumors)), one_hot(tumors, tumor_levels)])
    qd = orthonormal_basis(d)
    yres = y - qd @ (qd.T @ y)
    sse0 = np.sum(yres ** 2, axis=0)

    def f_for(lbl: np.ndarray) -> tuple[np.ndarray, np.ndarray, int, int]:
        a = one_hot(lbl, levels)
        ar = a - qd @ (qd.T @ a)
        qa = orthonormal_basis(ar)
        explained = np.sum((qa.T @ yres) ** 2, axis=0)
        sse1 = np.maximum(sse0 - explained, 0.0)
        df1 = qa.shape[1]
        df2 = len(lbl) - qd.shape[1] - df1
        with np.errstate(divide="ignore", invalid="ignore"):
            f = (explained / df1) / (sse1 / df2)
            r2 = explained / sse0
        return f, r2, df1, df2

    f_obs, r2, df1, df2 = f_for(labels)
    ge = np.zeros(y.shape[1], dtype=int)
    valid = np.zeros(y.shape[1], dtype=int)
    for _ in range(Q1_NPERM):
        f_perm, _, _, _ = f_for(permute_within_tumor(labels, tumors, rng))
        ok = np.isfinite(f_perm) & np.isfinite(f_obs)
        valid += ok
        ge += ok & (f_perm >= f_obs - 1e-12)
    p_emp = (ge + 1) / (valid + 1)
    return pd.DataFrame({
        "partial_r2": np.clip(r2, 0, 1),
        "omnibus_f": f_obs,
        "df_anatomy": df1,
        "df_residual": df2,
        "parametric_p": stats.f.sf(f_obs, df1, df2),
        "permutation_p": p_emp,
        "permutation_valid_n": valid,
        "bh_fdr": bh(p_emp),
    })

def grouped_predictions(y: np.ndarray, labels: np.ndarray, tumors: np.ndarray,
                        levels: list[str]) -> np.ndarray:
    n, m = y.shape
    pred = np.full((n, m), -1, dtype=np.int16)
    for tumor in pd.unique(tumors):
        test = tumors == tumor
        train = ~test
        centroids = np.full((len(levels), m), np.nan)
        for c, level in enumerate(levels):
            rows = train & (labels == level)
            if rows.any():
                centroids[c] = np.mean(y[rows], axis=0)
        distances = np.abs(y[test, None, :] - centroids[None, :, :])
        distances = np.where(np.isfinite(distances), distances, np.inf)
        pred[test] = np.argmin(distances, axis=1)
    return pred

def classification_metrics(pred: np.ndarray, labels: np.ndarray,
                           levels: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    true = np.array([levels.index(x) for x in labels], dtype=np.int16)
    m = pred.shape[1]
    recalls = np.full((len(levels), m), np.nan)
    f1s = np.full((len(levels), m), np.nan)
    confusion = np.zeros((len(levels), len(levels), m), dtype=int)
    for c in range(len(levels)):
        true_c = true == c
        for d in range(len(levels)):
            confusion[c, d] = np.sum(true_c[:, None] & (pred == d), axis=0)
        tp = confusion[c, c].astype(float)
        n_true = true_c.sum()
        n_pred = np.sum(pred == c, axis=0)
        recalls[c] = tp / n_true if n_true else np.nan
        denom = n_true + n_pred
        f1s[c] = np.divide(2 * tp, denom, out=np.zeros_like(tp), where=denom > 0)
    return np.nanmean(recalls, axis=0), np.nanmean(f1s, axis=0), confusion

if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("input");p.add_argument("--output",required=True)
    a=p.parse_args();z=np.load(a.input,allow_pickle=False)
    y=z["scores"];labels=z["labels"].astype(str);tumors=z["groups"].astype(str)
    if y.shape!=(122,201) or len(set(tumors))!=10 or not np.isfinite(y).all():raise ValueError("Frozen Ivy geometry mismatch")
    rng=np.random.default_rng(2026083101); result={}
    for branch,levels in BRANCHES.items():
        keep=np.isin(labels,levels);yy=y[keep];ll=labels[keep];tt=tumors[keep]
        assoc=q1_informativeness(yy,ll,tt,levels,rng)
        ba,_,_=classification_metrics(grouped_predictions(yy,ll,tt,levels),ll,levels)
        result[branch]={"median_partial_r2":float(assoc.partial_r2.median()),"associated_n":int((assoc.bh_fdr<0.05).sum()),"median_ba":float(np.median(ba)),"ba_ge_0_60":int((ba>=0.60).sum())}
    Path(a.output).write_text(json.dumps(result,indent=2)+"\n");print(json.dumps(result))
