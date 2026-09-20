#!/usr/bin/env python3
"""Prespecified R10A patient-level statistical analysis.

Requires the frozen pre-association endpoint receipt, frozen BayesPrism reconstruction,
and outcome-blind optimal MES matches. No publication figure is generated.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd


DATA_ROOT = Path(os.environ.get("GBM_MES_DATA_ROOT", "data"))
RESULT_ROOT = Path(os.environ.get("GBM_MES_OUTPUT_ROOT", "results"))
INPUT = DATA_ROOT / "prepared" / "cptac"
OUT = RESULT_ROOT / "validation" / "cptac"
SEED = 20260723
CV_REPEATS = 100
CV_FOLDS = 5
BOOTSTRAPS = 2000
PERMUTATIONS = 2000
MATCH_NULL = 5000


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(df: pd.DataFrame, name: str) -> None:
    df.to_csv(OUT / name, sep="\t", index=False, na_rep="NA")


def bh(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    x = p[ok]
    order = np.argsort(x)
    ranked = x[order]
    adj = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    tmp = np.empty_like(adj)
    tmp[order] = np.minimum(adj, 1.0)
    q[ok] = tmp
    return q


def parse_col(col: str) -> tuple[str, str, str]:
    f, estimator, endpoint = col.split("|", 2)
    return f, estimator, endpoint


def make_folds(meta: pd.DataFrame, repeats: int = CV_REPEATS) -> np.ndarray:
    rng = np.random.default_rng(SEED)
    folds = np.empty((repeats, len(meta)), dtype=np.int8)
    plex = meta["tmt_plex"].astype(str).to_numpy()
    for r in range(repeats):
        current = np.empty(len(meta), dtype=np.int8)
        for value in sorted(np.unique(plex)):
            idx = np.where(plex == value)[0]
            idx = rng.permutation(idx)
            offset = int(rng.integers(CV_FOLDS))
            current[idx] = (np.arange(len(idx)) + offset) % CV_FOLDS
        folds[r] = current
    return folds


def fixed_plex_dummies(meta: pd.DataFrame) -> np.ndarray:
    levels = sorted(meta["tmt_plex"].astype(str).unique())
    arr = meta["tmt_plex"].astype(str).to_numpy()
    return np.column_stack([(arr == level).astype(float) for level in levels[1:]])


def fold_design(meta: pd.DataFrame, train: np.ndarray, test: np.ndarray, extras: np.ndarray | None) -> tuple[np.ndarray, np.ndarray]:
    base = meta[["quality_score", "retained_missing_fraction"]].to_numpy(float)
    cont = base if extras is None or extras.shape[1] == 0 else np.column_stack([base, extras])
    mu = np.nanmean(cont[train], axis=0)
    sd = np.nanstd(cont[train], axis=0, ddof=1)
    sd[~np.isfinite(sd) | (sd == 0)] = 1.0
    tr = (np.where(np.isfinite(cont[train]), cont[train], mu) - mu) / sd
    te = (np.where(np.isfinite(cont[test]), cont[test], mu) - mu) / sd
    dummy = fixed_plex_dummies(meta)
    return np.column_stack([np.ones(len(train)), dummy[train], tr]), np.column_stack([np.ones(len(test)), dummy[test], te])


def cv_errors(y: np.ndarray, meta: pd.DataFrame, extras: np.ndarray | None, folds: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n, e = y.shape
    errors = np.full((folds.shape[0], n, e), np.nan, dtype=np.float32)
    baseline = np.full_like(errors, np.nan)
    predictions = np.full_like(errors, np.nan)
    for r in range(folds.shape[0]):
        for f in range(CV_FOLDS):
            test = np.where(folds[r] == f)[0]
            train = np.where(folds[r] != f)[0]
            xtr, xte = fold_design(meta, train, test, extras)
            ytr = y[train]
            mu = np.nanmean(ytr, axis=0)
            sd = np.nanstd(ytr, axis=0, ddof=1)
            valid_endpoint = np.isfinite(sd) & (sd > 0)
            sd[~valid_endpoint] = 1.0
            ytrz = (ytr - mu) / sd
            ytrz[~np.isfinite(ytrz)] = 0.0
            coef = np.linalg.lstsq(xtr, ytrz, rcond=None)[0]
            pred = xte @ coef
            truth = (y[test] - mu) / sd
            mask = np.isfinite(truth) & valid_endpoint[None, :]
            sqe = (truth - pred) ** 2
            sqb = truth ** 2
            errors[r, test] = np.where(mask, sqe, np.nan)
            baseline[r, test] = np.where(mask, sqb, np.nan)
            predictions[r, test] = np.where(mask, pred, np.nan)
    return errors, baseline, predictions


def q2_from_errors(errors: np.ndarray, baseline: np.ndarray) -> np.ndarray:
    sse = np.nansum(errors, axis=(0, 1))
    ssb = np.nansum(baseline, axis=(0, 1))
    return np.where(ssb > 0, 1.0 - sse / ssb, np.nan)


def q2_by_repeat(errors: np.ndarray, baseline: np.ndarray) -> np.ndarray:
    sse = np.nansum(errors, axis=1)
    ssb = np.nansum(baseline, axis=1)
    return np.where(ssb > 0, 1.0 - sse / ssb, np.nan)


def bootstrap_q2(errors: np.ndarray, baseline: np.ndarray, rng: np.random.Generator, b: int = BOOTSTRAPS) -> np.ndarray:
    n = errors.shape[1]
    out = np.empty((b, errors.shape[2]), dtype=np.float32)
    for i in range(b):
        idx = rng.integers(0, n, n)
        sse = np.nansum(errors[:, idx], axis=(0, 1))
        ssb = np.nansum(baseline[:, idx], axis=(0, 1))
        out[i] = np.where(ssb > 0, 1.0 - sse / ssb, np.nan)
    return out


def full_design(meta: pd.DataFrame, extras: np.ndarray | None = None) -> np.ndarray:
    train = np.arange(len(meta))
    x, _ = fold_design(meta, train, train, extras)
    return x


def residualize_vector(v: np.ndarray, x: np.ndarray) -> np.ndarray:
    ok = np.isfinite(v) & np.all(np.isfinite(x), axis=1)
    out = np.full(len(v), np.nan)
    if ok.sum() > x.shape[1] + 2:
        out[ok] = v[ok] - x[ok] @ np.linalg.lstsq(x[ok], v[ok], rcond=None)[0]
    return out


def partial_associations(y: np.ndarray, meta: pd.DataFrame, mes: np.ndarray, b: int = PERMUTATIONS) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(SEED + 11)
    x0 = full_design(meta, None)
    plex = meta["tmt_plex"].astype(str).to_numpy()
    betas = np.full(y.shape[1], np.nan)
    pvals = np.full(y.shape[1], np.nan)
    for j in range(y.shape[1]):
        ok = np.isfinite(y[:, j]) & np.isfinite(mes)
        if ok.sum() < x0.shape[1] + 5:
            continue
        ry = residualize_vector(y[:, j], x0)[ok]
        rm = residualize_vector(mes, x0)[ok]
        ry = (ry - ry.mean()) / ry.std(ddof=1)
        rm = (rm - rm.mean()) / rm.std(ddof=1)
        observed = float(np.mean(ry * rm))
        betas[j] = observed
        perm_stats = np.empty(b)
        plex_ok = plex[ok]
        for bi in range(b):
            rp = rm.copy()
            for level in np.unique(plex_ok):
                ii = np.where(plex_ok == level)[0]
                rp[ii] = rng.permutation(rp[ii])
            perm_stats[bi] = np.mean(ry * rp)
        pvals[j] = (1 + np.sum(np.abs(perm_stats) >= abs(observed))) / (b + 1)
    return betas, pvals


def endpoint_qvalues(pvals: np.ndarray, columns: list[str]) -> np.ndarray:
    q = np.full(len(columns), np.nan)
    groups: dict[tuple[str, str], list[int]] = {}
    for i, col in enumerate(columns):
        family, estimator, _ = parse_col(col)
        groups.setdefault((family, estimator), []).append(i)
    for idx in groups.values():
        q[idx] = bh(pvals[idx])
    return q


def residual_structure(y: np.ndarray, columns: list[str], meta: pd.DataFrame, mes: np.ndarray, b: int = PERMUTATIONS) -> pd.DataFrame:
    rng = np.random.default_rng(SEED + 21)
    plex = meta["tmt_plex"].astype(str).to_numpy()
    x1 = full_design(meta, mes[:, None])
    groups = {
        "JOINT_A_C": [i for i, c in enumerate(columns) if parse_col(c)[0] in {"A", "C"}],
        "FAMILY_A": [i for i, c in enumerate(columns) if parse_col(c)[0] == "A"],
        "FAMILY_C": [i for i, c in enumerate(columns) if parse_col(c)[0] == "C"],
    }
    rows = []
    for name, idx in groups.items():
        rmat = np.full((len(meta), len(idx)), np.nan)
        for k, j in enumerate(idx):
            rmat[:, k] = residualize_vector(y[:, j], x1)
            sd = np.nanstd(rmat[:, k], ddof=1)
            if np.isfinite(sd) and sd > 0:
                rmat[:, k] = (rmat[:, k] - np.nanmean(rmat[:, k])) / sd
        complete = np.where(np.isfinite(rmat), rmat, 0.0)
        eig = np.linalg.eigvalsh(np.cov(complete, rowvar=False))
        observed = float(eig[-1] / eig.sum())
        null = np.empty(b)
        for bi in range(b):
            z = rmat.copy()
            for k in range(z.shape[1]):
                for level in np.unique(plex):
                    ii = np.where((plex == level) & np.isfinite(z[:, k]))[0]
                    z[ii, k] = rng.permutation(z[ii, k])
            zz = np.where(np.isfinite(z), z, 0.0)
            ev = np.linalg.eigvalsh(np.cov(zz, rowvar=False))
            null[bi] = ev[-1] / ev.sum()
        rows.append({
            "record_type": "CONDITIONAL_RESIDUAL_STRUCTURE", "distance_family": name, "n_patients": len(meta),
            "n_endpoints": len(idx), "observed": observed, "null_median": float(np.median(null)),
            "effect": observed - float(np.median(null)), "permutation_p": (1 + np.sum(null >= observed)) / (b + 1),
            "permutations": b, "interpretation": "structured residual phospho-signaling; no rank claim",
        })
    out = pd.DataFrame(rows)
    out["BH_q"] = bh(out["permutation_p"].to_numpy())
    return out


def conditional_axis_permutation(y: np.ndarray, columns: list[str], meta: pd.DataFrame, mes: np.ndarray,
                                 axes: np.ndarray, folds: np.ndarray, m1_errors: np.ndarray,
                                 baseline: np.ndarray, observed_delta: np.ndarray, b: int = PERMUTATIONS) -> np.ndarray:
    rng = np.random.default_rng(SEED + 31)
    x1 = full_design(meta, mes[:, None])
    fitted = x1 @ np.linalg.lstsq(x1, axes, rcond=None)[0]
    residual = axes - fitted
    plex = meta["tmt_plex"].astype(str).to_numpy()
    exceed = np.zeros(y.shape[1], dtype=int)
    q2_m1 = q2_from_errors(m1_errors, baseline)
    for bi in range(b):
        rp = residual.copy()
        for level in np.unique(plex):
            ii = np.where(plex == level)[0]
            rp[ii] = rp[rng.permutation(ii)]
        perm_axes = fitted + rp
        e2, b2, _ = cv_errors(y, meta, np.column_stack([mes, perm_axes]), folds)
        delta = q2_from_errors(e2, b2) - q2_m1
        exceed += np.where(np.isfinite(delta) & np.isfinite(observed_delta) & (delta >= observed_delta), 1, 0)
        if (bi + 1) % 100 == 0:
            print(json.dumps({"axis_permutations_complete": bi + 1}), flush=True)
    return (1 + exceed) / (b + 1)


def standardize_activity(y: np.ndarray) -> np.ndarray:
    mu = np.nanmean(y, axis=0)
    sd = np.nanstd(y, axis=0, ddof=1)
    sd[~np.isfinite(sd) | (sd == 0)] = 1.0
    return (y - mu) / sd


def pair_distance(z: np.ndarray, i: np.ndarray, j: np.ndarray, min_fraction: float = 0.70) -> np.ndarray:
    a, b = z[i], z[j]
    ok = np.isfinite(a) & np.isfinite(b)
    n = ok.sum(axis=1)
    d2 = np.where(ok, (a - b) ** 2, 0.0).sum(axis=1)
    return np.where(n >= math.ceil(min_fraction * z.shape[1]), np.sqrt(d2 / np.maximum(n, 1)), np.nan)


def constrained_pair_null(pairs: pd.DataFrame, meta: pd.DataFrame, b: int = MATCH_NULL) -> tuple[np.ndarray, float]:
    rng = np.random.default_rng(SEED + 41)
    id_to_i = {x: i for i, x in enumerate(meta["case_id"])}
    left = np.array([id_to_i[x] for x in pairs["case_i"]], dtype=int)
    right = np.array([id_to_i[x] for x in pairs["case_j"]], dtype=int)
    mes = meta["bulk_mes_z"].to_numpy(float)
    quality = meta["quality_score"].to_numpy(float)
    plex = meta["tmt_plex"].astype(str).to_numpy()
    qdiff = []
    for i in range(len(meta)):
        for j in range(i + 1, len(meta)):
            if abs(mes[i] - mes[j]) <= 0.20:
                qdiff.append(abs(quality[i] - quality[j]))
    qcuts = np.quantile(qdiff, [0.25, 0.50, 0.75])

    def slot_group(i: int, j: int) -> tuple[int, int, int]:
        dm = min(3, int(abs(mes[i] - mes[j]) / 0.05))
        dq = int(np.searchsorted(qcuts, abs(quality[i] - quality[j]), side="right"))
        return int(plex[i] != plex[j]), dm, dq

    targets = [slot_group(i, j) for i, j in zip(left, right)]
    current = right.copy()
    configs = np.empty((b, len(current)), dtype=np.int16)
    accepted = 0
    attempted = 0
    warmup = 2000
    total_samples = b + 1
    for step in range(warmup + total_samples * 100):
        a, c = rng.choice(len(current), 2, replace=False)
        ja, jc = current[a], current[c]
        attempted += 1
        valid = left[a] != jc and left[c] != ja
        valid = valid and abs(mes[left[a]] - mes[jc]) <= 0.20 and abs(mes[left[c]] - mes[ja]) <= 0.20
        valid = valid and slot_group(left[a], jc) == targets[a] and slot_group(left[c], ja) == targets[c]
        if valid:
            current[a], current[c] = jc, ja
            accepted += 1
        # Sample the chain on the frozen attempt schedule even when the current
        # proposal is rejected.  Sampling only on accepted proposals at exact
        # schedule steps left uninitialized rows in the configuration matrix.
        if step >= warmup and (step - warmup) % 100 == 0:
            bi = (step - warmup) // 100
            if bi < b:
                configs[bi] = current
    if np.any(configs < 0) or np.any(configs >= len(meta)):
        raise RuntimeError("Constrained matched-null configuration contains an invalid patient index")
    return np.column_stack([np.tile(left, (b, 1)), configs]), accepted / max(attempted, 1)


def matched_summary(y: np.ndarray, columns: list[str], pairs: pd.DataFrame, meta: pd.DataFrame,
                    null_configs: np.ndarray, sensitivity: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    id_to_i = {x: i for i, x in enumerate(meta["case_id"])}
    i = np.array([id_to_i[x] for x in pairs["case_i"]], dtype=int)
    j = np.array([id_to_i[x] for x in pairs["case_j"]], dtype=int)
    groups = {
        "JOINT_A_C": [k for k, c in enumerate(columns) if parse_col(c)[0] in {"A", "C"}],
        "FAMILY_A": [k for k, c in enumerate(columns) if parse_col(c)[0] == "A"],
        "FAMILY_C": [k for k, c in enumerate(columns) if parse_col(c)[0] == "C"],
        "UNBIASED_KSEA": [k for k, c in enumerate(columns) if parse_col(c)[0] == "B" and parse_col(c)[1] == "KSEA_OMNIPATH_CURATED_ZMEAN"],
    }
    rows, pair_rows = [], []
    for name, idx in groups.items():
        z = standardize_activity(y[:, idx])
        observed_pair = pair_distance(z, i, j)
        observed = float(np.nanmedian(observed_pair))
        null = np.empty(len(null_configs))
        m = len(i)
        for b in range(len(null_configs)):
            null[b] = np.nanmedian(pair_distance(z, null_configs[b, :m], null_configs[b, m:]))
        rows.append({
            "record_type": "SUMMARY", "sensitivity": sensitivity, "distance_family": name,
            "n_pairs": int(np.isfinite(observed_pair).sum()), "n_endpoints": len(idx), "observed_median_distance": observed,
            "null_median": float(np.median(null)), "effect": observed - float(np.median(null)),
            "ratio": observed / float(np.median(null)) if np.median(null) > 0 else np.nan,
            "permutation_p": (1 + np.sum(null >= observed)) / (len(null) + 1), "BH_q": np.nan,
            "null_reassignments": len(null), "minimum_pairs_gate": 20, "status": "EVALUATED",
        })
        if name == "JOINT_A_C" and sensitivity == "PRIMARY":
            for k, d in enumerate(observed_pair):
                pair_rows.append({
                    "record_type": "PAIR", "sensitivity": sensitivity, "distance_family": name,
                    "pair_id": pairs.iloc[k]["pair_id"], "case_i": pairs.iloc[k]["case_i"], "case_j": pairs.iloc[k]["case_j"],
                    "delta_mes": pairs.iloc[k]["delta_mes"], "plex_mismatch": pairs.iloc[k]["plex_mismatch"],
                    "n_pairs": 1, "n_endpoints": len(idx), "observed_median_distance": d, "null_median": np.nan,
                    "effect": np.nan, "ratio": np.nan, "permutation_p": np.nan, "BH_q": np.nan,
                    "null_reassignments": len(null), "minimum_pairs_gate": 20, "status": "EVALUATED_PAIR",
                })
    summary = pd.DataFrame(rows)
    summary["BH_q"] = bh(summary["permutation_p"].to_numpy())
    return summary, pd.DataFrame(pair_rows)


def sensitivity_cv(y: np.ndarray, columns: list[str], meta: pd.DataFrame, mes: np.ndarray, axes: np.ndarray,
                   folds: np.ndarray, rng: np.random.Generator, label: str) -> pd.DataFrame:
    idx = [i for i, c in enumerate(columns) if parse_col(c)[0] in {"A", "C"}]
    yy = y[:, idx]
    e1, b1, _ = cv_errors(yy, meta, mes[:, None], folds)
    e2, b2, _ = cv_errors(yy, meta, np.column_stack([mes, axes]), folds)
    q1, q2 = q2_from_errors(e1, b1), q2_from_errors(e2, b2)
    boot1 = bootstrap_q2(e1, b1, rng)
    boot2 = bootstrap_q2(e2, b2, rng)
    delta_boot = boot2 - boot1
    beta, p = partial_associations(yy, meta, mes)
    q = endpoint_qvalues(p, [columns[i] for i in idx])
    rows = []
    for k, j in enumerate(idx):
        family, estimator, endpoint = parse_col(columns[j])
        rows.append({
            "sensitivity": label, "family": family, "estimator": estimator, "endpoint": endpoint,
            "n": int(np.isfinite(yy[:, k]).sum()), "MES_partial_beta": beta[k], "MES_permutation_p": p[k], "MES_BH_q": q[k],
            "Q2_M1": q1[k], "Q2_M2": q2[k], "delta_Q2": q2[k] - q1[k],
            "delta_Q2_boot_CI_low": np.nanquantile(delta_boot[:, k], 0.025),
            "delta_Q2_boot_CI_high": np.nanquantile(delta_boot[:, k], 0.975),
            "positive_delta_repeat_fraction": float(np.mean((q2_by_repeat(e2, b2)[:, k] - q2_by_repeat(e1, b1)[:, k]) > 0)),
        })
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    endpoint_hash = sha256(INPUT / "R10A_SIGNALING_ENDPOINT_FREEZE.tsv")
    if endpoint_hash != "765e3277cf09e48df1aa269aa30906c87d7f4193e90ffd20fe2161cf2a271285":
        raise RuntimeError("Endpoint freeze hash mismatch")
    pre = json.loads((INPUT / "R10A_PREASSOCIATION_RECEIPT.json").read_text())
    bp = json.loads((INPUT / "R10A_BAYESPRISM_RECEIPT.json").read_text())
    if pre["association_outcomes_inspected"] or bp["signaling_outcomes_inspected"]:
        raise RuntimeError("Pre-association freeze invalid")

    meta = pd.read_csv(INPUT / "R10A_PATIENT_INPUT_MATRIX.tsv", sep="\t")
    activity = pd.read_csv(INPUT / "R10A_SIGNALING_ACTIVITY_MATRIX.tsv", sep="\t").set_index("case_id").loc[meta["case_id"]]
    adjusted = pd.read_csv(INPUT / "R10A_SIGNALING_ACTIVITY_TOTAL_PROTEIN_ADJUSTED.tsv", sep="\t").set_index("case_id").loc[meta["case_id"]]
    tautology = pd.read_csv(INPUT / "R10A_SIGNALING_ACTIVITY_TAUTOLOGY_EXCLUDED.tsv", sep="\t").set_index("case_id").loc[meta["case_id"]]
    columns = activity.columns.tolist()
    y = activity.to_numpy(float)
    y_adj = adjusted.to_numpy(float)
    y_taut = tautology.to_numpy(float)
    mes = meta["bulk_mes_z"].to_numpy(float)
    axes = meta[["malignant_mes_z", "myeloid_z", "hypoxia_pruned_z", "matrix_pruned_z", "vascular_stromal_z"]].to_numpy(float)
    if len(meta) != 98 or not np.isfinite(axes).all():
        raise RuntimeError("R10A final patient/realization matrix invalid")
    pairs = pd.read_csv(INPUT / "R10A_OPTIMAL_MES_MATCHED_PAIRS.tsv", sep="\t")
    if len(pairs) < 20:
        raise RuntimeError("Matched-pair gate did not pass")

    rng = np.random.default_rng(SEED + 51)
    folds = make_folds(meta)
    fold_rows = []
    for r in range(CV_REPEATS):
        for f in range(CV_FOLDS):
            for i in np.where(folds[r] == f)[0]:
                fold_rows.append((r + 1, f + 1, meta.iloc[i]["case_id"], meta.iloc[i]["tmt_plex"]))
    write_tsv(pd.DataFrame(fold_rows, columns=["repeat", "fold", "case_id", "tmt_plex"]), "R10A_REPEATED_CV_FOLDS.tsv")

    e0, b0, _ = cv_errors(y, meta, None, folds)
    e1, b1, p1 = cv_errors(y, meta, mes[:, None], folds)
    e2, b2, p2 = cv_errors(y, meta, np.column_stack([mes, axes]), folds)
    q0, q1, q2 = q2_from_errors(e0, b0), q2_from_errors(e1, b1), q2_from_errors(e2, b2)
    q1r, q2r = q2_by_repeat(e1, b1), q2_by_repeat(e2, b2)
    boot1 = bootstrap_q2(e1, b1, rng)
    boot2 = bootstrap_q2(e2, b2, rng)
    delta_boot = boot2 - boot1
    beta, assoc_p = partial_associations(y, meta, mes)
    assoc_q = endpoint_qvalues(assoc_p, columns)

    pred_rows = []
    for i, col in enumerate(columns):
        family, estimator, endpoint = parse_col(col)
        pred_rows.append({
            "record_type": "ENDPOINT", "family": family, "estimator": estimator, "endpoint": endpoint,
            "n": int(np.isfinite(y[:, i]).sum()), "MES_partial_beta": beta[i], "MES_permutation_p": assoc_p[i], "MES_BH_q": assoc_q[i],
            "Q2_M0": q0[i], "Q2_M1": q1[i], "delta_Q2_M1_vs_M0": q1[i] - q0[i],
            "Q2_M1_boot_CI_low": np.nanquantile(boot1[:, i], .025), "Q2_M1_boot_CI_high": np.nanquantile(boot1[:, i], .975),
            "RMSE_M1_standardized": float(np.sqrt(np.nanmean(e1[:, :, i]))),
            "informative_gate": "SUPPORTED" if assoc_q[i] < .05 or np.nanquantile(boot1[:, i], .025) > 0 else "NOT_SUPPORTED",
        })
    predictability = pd.DataFrame(pred_rows)
    for family in ["A", "B", "C"]:
        ii = np.array([i for i, c in enumerate(columns) if parse_col(c)[0] == family])
        if len(ii):
            predictability.loc[len(predictability)] = {
                "record_type": "FAMILY_MEDIAN", "family": family, "estimator": "ALL_RETAINED", "endpoint": "FAMILY_MEDIAN",
                "n": len(meta), "Q2_M0": np.nanmedian(q0[ii]), "Q2_M1": np.nanmedian(q1[ii]),
                "delta_Q2_M1_vs_M0": np.nanmedian(q1[ii] - q0[ii]),
                "informative_gate": "DESCRIPTIVE_FAMILY_SUMMARY",
            }
    write_tsv(predictability, "R10A_MES_SIGNALING_PREDICTABILITY.tsv")

    observed_delta = q2 - q1
    perm_p = conditional_axis_permutation(y, columns, meta, mes, axes, folds, e1, b1, observed_delta)
    perm_q = endpoint_qvalues(perm_p, columns)
    incr_rows = []
    for i, col in enumerate(columns):
        family, estimator, endpoint = parse_col(col)
        ci_low = np.nanquantile(delta_boot[:, i], .025)
        ci_high = np.nanquantile(delta_boot[:, i], .975)
        pos = float(np.mean((q2r[:, i] - q1r[:, i]) > 0))
        supported = observed_delta[i] > 0 and ci_low > 0 and perm_q[i] < .05 and pos >= .80
        incr_rows.append({
            "record_type": "ENDPOINT", "family": family, "estimator": estimator, "endpoint": endpoint,
            "n": int(np.isfinite(y[:, i]).sum()), "Q2_M1": q1[i], "Q2_M2": q2[i], "delta_Q2": observed_delta[i],
            "delta_Q2_boot_CI_low": ci_low, "delta_Q2_boot_CI_high": ci_high,
            "conditional_permutation_p": perm_p[i], "BH_q": perm_q[i], "positive_delta_repeat_fraction": pos,
            "incremental_gate": "SUPPORTED" if supported else "NOT_SUPPORTED",
            "realization_axes": "malignant_mes_z;myeloid_z;hypoxia_pruned_z;matrix_pruned_z;vascular_stromal_z",
        })
    incremental = pd.DataFrame(incr_rows)
    for family in ["A", "B", "C"]:
        ii = np.array([i for i, c in enumerate(columns) if parse_col(c)[0] == family])
        if len(ii):
            incremental.loc[len(incremental)] = {
                "record_type": "FAMILY_MEDIAN", "family": family, "estimator": "ALL_RETAINED", "endpoint": "FAMILY_MEDIAN",
                "n": len(meta), "Q2_M1": np.nanmedian(q1[ii]), "Q2_M2": np.nanmedian(q2[ii]),
                "delta_Q2": np.nanmedian(observed_delta[ii]), "incremental_gate": "DESCRIPTIVE_FAMILY_SUMMARY",
                "realization_axes": "five frozen axes",
            }
    write_tsv(incremental, "R10A_INCREMENTAL_REALIZATION_VALUE.tsv")

    structure = residual_structure(y, columns, meta, mes)
    null_configs, acceptance = constrained_pair_null(pairs, meta)
    primary_match, pair_rows = matched_summary(y, columns, pairs, meta, null_configs, "PRIMARY")
    primary_match["null_chain_swap_acceptance"] = acceptance
    matched_all = pd.concat([primary_match, pair_rows], ignore_index=True, sort=False)
    structure_rows = structure.copy()
    for col in matched_all.columns:
        if col not in structure_rows:
            structure_rows[col] = np.nan
    for col in structure_rows.columns:
        if col not in matched_all:
            matched_all[col] = np.nan
    matched_all = pd.concat([matched_all[structure_rows.columns], structure_rows], ignore_index=True)
    write_tsv(matched_all, "R10A_MATCHED_NON_EQUIVALENCE.tsv")

    adj_cv = sensitivity_cv(y_adj, columns, meta, mes, axes, folds, rng, "TOTAL_PROTEIN_ADJUSTED")
    adj_match, _ = matched_summary(y_adj, columns, pairs, meta, null_configs, "TOTAL_PROTEIN_ADJUSTED")
    adj_match["result_type"] = "MATCHED_NON_EQUIVALENCE"
    adj_cv["result_type"] = "ENDPOINT_CV_AND_ASSOCIATION"
    for col in set(adj_match.columns) - set(adj_cv.columns): adj_cv[col] = np.nan
    for col in set(adj_cv.columns) - set(adj_match.columns): adj_match[col] = np.nan
    protein_results = pd.concat([adj_cv, adj_match[adj_cv.columns]], ignore_index=True)
    write_tsv(protein_results, "R10A_TOTAL_PROTEIN_ADJUSTED_RESULTS.tsv")

    taut_cv = sensitivity_cv(y_taut, columns, meta, mes, axes, folds, rng, "CANONICAL_MES_GENE_PROTEINS_EXCLUDED")
    taut_match, _ = matched_summary(y_taut, columns, pairs, meta, null_configs, "CANONICAL_MES_GENE_PROTEINS_EXCLUDED")
    membership = pd.read_csv(INPUT / "R10A_ENDPOINT_SITE_MEMBERSHIP.tsv", sep="\t")
    mes_map = pd.read_csv(INPUT / "R10A_MES_MAPPING_AUDIT.tsv", sep="\t")
    mes_genes = set(mes_map.loc[mes_map["record_type"].eq("GENE"), "gene"])
    membership["gene"] = membership["site_id"].str.rsplit("_", n=1).str[0]
    overlap = membership.loc[membership["mapped_primary_site"].astype(str).str.upper().eq("TRUE")].groupby(["family", "estimator", "endpoint"])["gene"].apply(lambda x: len(set(x) & mes_genes)).reset_index(name="mapped_canonical_MES_proteins")
    taut_cv = taut_cv.merge(overlap, on=["family", "estimator", "endpoint"], how="left")
    taut_cv["result_type"] = "ENDPOINT_CV_AND_ASSOCIATION"
    taut_match["result_type"] = "MATCHED_NON_EQUIVALENCE"
    for col in set(taut_match.columns) - set(taut_cv.columns): taut_cv[col] = np.nan
    for col in set(taut_cv.columns) - set(taut_match.columns): taut_match[col] = np.nan
    taut_results = pd.concat([taut_cv, taut_match[taut_cv.columns]], ignore_index=True)
    write_tsv(taut_results, "R10A_TAUTOLOGY_AUDIT.tsv")

    # Cross-estimator concordance for kinase names shared by primary and secondary resources.
    def canonical_kinase(name: str) -> str:
        return name.split("/")[-1].split(".")[0].upper()
    ksea = predictability[(predictability["record_type"] == "ENDPOINT") & (predictability["family"] == "B") & (predictability["estimator"] == "KSEA_OMNIPATH_CURATED_ZMEAN")].copy()
    ptm = predictability[(predictability["record_type"] == "ENDPOINT") & (predictability["family"] == "B") & (predictability["estimator"] == "PTMSIGDB_SIGNED_ZMEAN")].copy()
    ksea["canonical_kinase"] = ksea["endpoint"].map(canonical_kinase)
    ptm["canonical_kinase"] = ptm["endpoint"].map(canonical_kinase)
    concord = ksea.merge(ptm, on="canonical_kinase", suffixes=("_KSEA", "_PTMSIGDB"))
    concord["MES_effect_direction_concordant"] = np.sign(concord["MES_partial_beta_KSEA"]) == np.sign(concord["MES_partial_beta_PTMSIGDB"])
    concord["both_BH_q_lt_0_05"] = (concord["MES_BH_q_KSEA"] < .05) & (concord["MES_BH_q_PTMSIGDB"] < .05)
    write_tsv(concord, "R10A_ACTIVITY_ESTIMATOR_CONCORDANCE.tsv")

    # Influence audit targets prespecified primary matched distance and the strongest formal incremental endpoint.
    primary_joint = primary_match.loc[primary_match["distance_family"].eq("JOINT_A_C")].iloc[0]
    z_joint = standardize_activity(y[:, [i for i, c in enumerate(columns) if parse_col(c)[0] in {"A", "C"}]])
    id_to_i = {x: i for i, x in enumerate(meta["case_id"])}
    pi = np.array([id_to_i[x] for x in pairs["case_i"]]); pj = np.array([id_to_i[x] for x in pairs["case_j"]])
    pdist = pair_distance(z_joint, pi, pj)
    loo_match_effect = []
    for case in meta["case_id"]:
        keep = ~(pairs["case_i"].eq(case) | pairs["case_j"].eq(case)).to_numpy()
        loo_match_effect.append(np.nanmedian(pdist[keep]) - primary_joint["null_median"])
    formal = incremental[(incremental["record_type"] == "ENDPOINT") & incremental["family"].isin(["A", "C"])].copy()
    best = formal.sort_values("delta_Q2", ascending=False).iloc[0]
    best_col = columns.index(f"{best['family']}|{best['estimator']}|{best['endpoint']}")
    loo_delta = []
    for omit in range(len(meta)):
        keep = np.arange(len(meta)) != omit
        mm = meta.loc[keep].reset_index(drop=True)
        ff = make_folds(mm)
        yy = y[keep, best_col:best_col+1]
        ee1, bb1, _ = cv_errors(yy, mm, mes[keep, None], ff)
        ee2, bb2, _ = cv_errors(yy, mm, np.column_stack([mes[keep], axes[keep]]), ff)
        loo_delta.append(float(q2_from_errors(ee2, bb2)[0] - q2_from_errors(ee1, bb1)[0]))
    influence = pd.DataFrame([
        {"target": "JOINT_A_C_MATCHED_EFFECT", "endpoint": "JOINT_A_C", "full_estimate": primary_joint["effect"], "LOO_min": np.nanmin(loo_match_effect), "LOO_max": np.nanmax(loo_match_effect), "sign_flip": bool(np.nanmin(loo_match_effect) <= 0), "status": "ROBUST" if np.nanmin(loo_match_effect) > 0 else "OUTLIER_DEPENDENT"},
        {"target": "STRONGEST_INCREMENTAL_DELTA_Q2", "endpoint": best["endpoint"], "full_estimate": best["delta_Q2"], "LOO_min": np.nanmin(loo_delta), "LOO_max": np.nanmax(loo_delta), "sign_flip": bool(np.nanmin(loo_delta) <= 0), "status": "ROBUST" if np.nanmin(loo_delta) > 0 else "OUTLIER_DEPENDENT"},
    ])
    write_tsv(influence, "R10A_INFLUENCE_AUDIT.tsv")

    matched_supported = bool(((primary_match["distance_family"].isin(["JOINT_A_C", "FAMILY_A", "FAMILY_C"])) & (primary_match["BH_q"] < .05)).any())
    residual_supported = bool((structure.loc[structure["distance_family"].eq("JOINT_A_C"), "BH_q"] < .05).any())
    protein_matched_supported = bool(((adj_match["distance_family"].isin(["JOINT_A_C", "FAMILY_A", "FAMILY_C"])) & (adj_match["BH_q"] < .05)).any())
    taut_matched_supported = bool(((taut_match["distance_family"].isin(["JOINT_A_C", "FAMILY_A", "FAMILY_C"])) & (taut_match["BH_q"] < .05)).any())
    influence_supported = bool((influence["status"] == "ROBUST").all())
    nonident = matched_supported and residual_supported and protein_matched_supported and taut_matched_supported and influence_supported

    supported_endpoints = incremental[(incremental["record_type"] == "ENDPOINT") & incremental["family"].isin(["A", "C"]) & incremental["incremental_gate"].eq("SUPPORTED")]
    incremental_supported = False
    robust_incremental_endpoints = []
    for _, row in supported_endpoints.iterrows():
        ap = adj_cv.loc[(adj_cv["family"] == row["family"]) & (adj_cv["endpoint"] == row["endpoint"])]
        tp = taut_cv.loc[(taut_cv["family"] == row["family"]) & (taut_cv["endpoint"] == row["endpoint"])]
        if len(ap) and len(tp) and ap.iloc[0]["delta_Q2_boot_CI_low"] > 0 and tp.iloc[0]["delta_Q2_boot_CI_low"] > 0:
            incremental_supported = True
            robust_incremental_endpoints.append(row["endpoint"])

    mes_informative = bool(((predictability["record_type"] == "ENDPOINT") & predictability["family"].isin(["A", "C"]) & (predictability["informative_gate"] == "SUPPORTED")).any())
    if nonident and incremental_supported:
        status = "FUNCTIONAL_NONIDENTIFICATION_WITH_INCREMENTAL_REALIZATION_SUPPORTED"
    elif nonident:
        status = "FUNCTIONAL_NONIDENTIFICATION_SUPPORTED"
    elif mes_informative:
        status = "MES_SIGNALING_ASSOCIATION_WITHOUT_NONIDENTIFICATION_UPGRADE"
    else:
        status = "PHOSPHOPROTEOMIC_NON_EQUIVALENCE_NOT_SUPPORTED"

    best_incr = incremental[incremental["record_type"].eq("ENDPOINT")].sort_values("delta_Q2", ascending=False).iloc[0]
    anchored_sig = predictability[(predictability["record_type"] == "ENDPOINT") & (predictability["family"] == "A") & (predictability["MES_BH_q"] < .05)]
    unbiased_sig = predictability[(predictability["record_type"] == "ENDPOINT") & (predictability["family"] == "B") & (predictability["estimator"] == "KSEA_OMNIPATH_CURATED_ZMEAN") & (predictability["MES_BH_q"] < .05)]
    interim = f"""# R10A interim adjudication

## Verdict

`{status}`

## Authenticated denominator and gates

- Phosphoproteomic source: `PDC000205`; `PDC000204` is the proteome study.
- Official tumor design: 100 patients; RNA available: 99; PDC000205 phosphosite available: 99; strict authenticated intersection and final evaluable denominator: **98 patients**.
- Canonical MES mapping: **95/95**, passing the frozen R10 gate of 86/95.
- Eligible primary phosphosites: **7,683** after single-site/unambiguous-gene, 70% coverage, and variance gates.
- MES-matched analysis: **{len(pairs)} unique non-overlapping pairs** at the frozen 0.20-SD caliper.

## Formal evidence

- MES-conditioned joint residual structure: observed PC1 variance fraction {structure.loc[structure.distance_family.eq('JOINT_A_C'),'observed'].iloc[0]:.4f}, BH q={structure.loc[structure.distance_family.eq('JOINT_A_C'),'BH_q'].iloc[0]:.4g}.
- Joint Family A+C matched distance: observed {primary_joint['observed_median_distance']:.4f} versus constrained-null median {primary_joint['null_median']:.4f}, effect {primary_joint['effect']:.4f}, BH q={primary_joint['BH_q']:.4g}.
- Total-protein-adjusted matched gate: {'supported' if protein_matched_supported else 'not supported'}; canonical-MES-protein-excluded gate: {'supported' if taut_matched_supported else 'not supported'}.
- Prespecified Family A MES associations surviving BH: {len(anchored_sig)}. Unbiased KSEA MES associations surviving BH: {len(unbiased_sig)}.
- Strongest incremental endpoint: `{best_incr['endpoint']}`, delta-Q2={best_incr['delta_Q2']:.4f} (95% patient-bootstrap CI {best_incr['delta_Q2_boot_CI_low']:.4f} to {best_incr['delta_Q2_boot_CI_high']:.4f}; conditional BH q={best_incr['BH_q']:.4g}).
- Incremental endpoints robust to total-protein and tautology sensitivities: {('; '.join(robust_incremental_endpoints)) if robust_incremental_endpoints else 'none'}.
- Influence audit: {', '.join(influence['status'])}.

## Interpretation boundary

This gate concerns patient-level cross-sectional phospho-signaling identity. It is not a kinase-causality claim, does not equate phosphosite abundance with activity without substrate inference, does not imply universal rank 2, and does not establish longitudinal rewiring. R10B proceeds independently regardless of this verdict; manuscript wording remains unchanged until cross-modal R10C.
"""
    (OUT / "R10A_INTERIM_ADJUDICATION.md").write_text(interim)

    receipt = {
        "status": status, "n_final": len(meta), "n_pairs": len(pairs), "matched_supported": matched_supported,
        "residual_structure_supported": residual_supported, "protein_adjusted_supported": protein_matched_supported,
        "tautology_excluded_supported": taut_matched_supported, "influence_supported": influence_supported,
        "incremental_supported": incremental_supported, "robust_incremental_endpoints": robust_incremental_endpoints,
        "strongest_incremental_endpoint": best_incr["endpoint"], "strongest_incremental_delta_Q2": float(best_incr["delta_Q2"]),
        "strongest_incremental_CI": [float(best_incr["delta_Q2_boot_CI_low"]), float(best_incr["delta_Q2_boot_CI_high"])],
        "strongest_incremental_BH_q": float(best_incr["BH_q"]), "endpoint_freeze_sha256": endpoint_hash,
    }
    (OUT / "R10A_FINAL_RECEIPT.json").write_text(json.dumps(receipt, indent=2) + "\n")
    np.savez_compressed(OUT / "R10A_OOF_NUMERICAL_RECEIPT.npz", folds=folds, errors_M0=e0, errors_M1=e1, errors_M2=e2,
                        baseline=b1, predictions_M1=p1, predictions_M2=p2, columns=np.array(columns))
    print(json.dumps(receipt, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
