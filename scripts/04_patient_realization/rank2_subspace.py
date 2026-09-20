#!/usr/bin/env python3
"""Formal R9A structural-validity and R9B measurement-invariance analyses."""

from __future__ import annotations

import json
import math
import os
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy


RESULT_ROOT = Path(os.environ.get("GBM_MES_OUTPUT_ROOT", "results"))
RUN = RESULT_ROOT / "construct_validation"
PATIENT_REALIZATION = RESULT_ROOT / "patient_realization"
AOUT = RUN / "01_structural_validity"
BOUT = RUN / "02_measurement_invariance"
LOG = RUN / "logs"
AXES = ["MALIGNANT_STATE", "MYELOID", "HYPOXIA", "MATRIX", "VASCULAR_STROMAL"]
RAW_COLS = ["malignant_mes_z", "myeloid_z", "hypoxia_z", "matrix_z", "vascular_z"]
RESID_COLS = [f"residual_z_{x}" for x in AXES]
COHORTS = ["CGGA325", "CGGA693"]


def write(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="NA", float_format="%.12g")


def eig(cov: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values, vectors = np.linalg.eigh(cov)
    order = np.argsort(values)[::-1]
    return values[order], vectors[:, order]


def covariance(x: np.ndarray) -> np.ndarray:
    return np.cov(x, rowvar=False, ddof=1)


def projector(u: np.ndarray, k: int) -> np.ndarray:
    q = u[:, :k]
    return q @ q.T


def projection_similarity(u: np.ndarray, v: np.ndarray, k: int) -> float:
    return float(np.trace(projector(u, k) @ projector(v, k)) / k)


def angles_deg(u: np.ndarray, v: np.ndarray, k: int) -> np.ndarray:
    s = np.linalg.svd(u[:, :k].T @ v[:, :k], compute_uv=False)
    return np.degrees(np.arccos(np.clip(s, -1.0, 1.0)))


def effective_rank(values: np.ndarray, p95a: np.ndarray, p95b: np.ndarray) -> int:
    rank = 0
    for j, value in enumerate(values):
        if value > p95a[j] and value > p95b[j]:
            rank += 1
        else:
            break
    return rank


def fit_basis(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    center = x.mean(axis=0)
    values, vectors = eig(covariance(x - center))
    return center, values, vectors


def reconstruct(x: np.ndarray, center: np.ndarray, vectors: np.ndarray, k: int) -> np.ndarray:
    z = x - center
    q = vectors[:, :k]
    return center + z @ q @ q.T


def q2_rmse(y: np.ndarray, pred: np.ndarray, baseline: np.ndarray | None = None) -> tuple[float, float]:
    if baseline is None:
        baseline = np.broadcast_to(y.mean(axis=0), y.shape)
    elif baseline.ndim == 1:
        baseline = np.broadcast_to(baseline, y.shape)
    sse = float(np.square(y - pred).sum())
    sst = float(np.square(y - baseline).sum())
    return 1.0 - sse / sst, math.sqrt(sse / y.size)


def folds(n: int, rng: np.random.Generator, k: int = 5) -> list[np.ndarray]:
    order = rng.permutation(n)
    return [a for a in np.array_split(order, k) if len(a)]


def slope_vector(s: np.ndarray, r: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.column_stack([np.ones(len(s)), s])
    coef = np.linalg.lstsq(x, r, rcond=None)[0]
    return coef[0], coef[1]


def residualize(s: np.ndarray, r: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    alpha, beta = slope_vector(s, r)
    return r - alpha - np.outer(s, beta), alpha, beta


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / denom) if denom > 0 else np.nan


def meaning_metrics(a: np.ndarray, b: np.ndarray) -> dict[str, float]:
    return {
        "cosine_similarity": cosine(a, b),
        "tucker_congruence": cosine(a, b),
        "standardized_euclidean_distance": float(np.sqrt(np.mean(np.square(a - b)))),
        "euclidean_distance": float(np.linalg.norm(a - b)),
    }


def rv_coefficient(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.trace(a @ b) / math.sqrt(np.trace(a @ a) * np.trace(b @ b)))


def matrix_correlation(a: np.ndarray, b: np.ndarray) -> float:
    ix = np.triu_indices_from(a)
    return float(np.corrcoef(a[ix], b[ix])[0, 1])


def geometry_metrics(a: np.ndarray, b: np.ndarray) -> dict[str, float]:
    va, ua = eig(a)
    vb, ub = eig(b)
    aa = angles_deg(ua, ub, 2)
    an = a / np.trace(a)
    bn = b / np.trace(b)
    return {
        "rv_coefficient": rv_coefficient(a, b),
        "matrix_correlation_upper": matrix_correlation(a, b),
        "trace_normalized_frobenius_distance": float(np.linalg.norm(an - bn, ord="fro")),
        "rank2_projection_similarity": projection_similarity(ua, ub, 2),
        "rank2_angle1_deg": float(aa[0]),
        "rank2_angle2_deg": float(aa[1]),
        "rank2_grassmann_distance_rad": float(np.linalg.norm(np.radians(aa))),
        "lambda1_ratio_325_to_693": float(va[0] / vb[0]),
    }


def pa_analysis(resid: dict[str, np.ndarray]) -> tuple[pd.DataFrame, dict[str, dict[str, np.ndarray]]]:
    rows: list[dict] = []
    thresholds: dict[str, dict[str, np.ndarray]] = {}
    for cohort in COHORTS:
        x = resid[cohort]
        n, p = x.shape
        obs, _ = eig(covariance(x))
        nulls: dict[str, np.ndarray] = {}
        for name, seed in [("INDEPENDENT_AXIS_PERMUTATION", 910101), ("GAUSSIAN_REFERENCE", 910102)]:
            rng = np.random.default_rng(seed + (0 if cohort == "CGGA325" else 10000))
            vals = np.empty((10000, p))
            variances = x.var(axis=0, ddof=1)
            for b in range(10000):
                if name == "INDEPENDENT_AXIS_PERMUTATION":
                    z = np.column_stack([x[rng.permutation(n), j] for j in range(p)])
                else:
                    z = rng.normal(size=(n, p)) * np.sqrt(variances)
                vals[b] = eig(covariance(z))[0]
            nulls[name] = vals
            for j in range(p):
                rows.append({
                    "cohort": cohort, "null_model": name, "component": j + 1, "n": n, "p": p,
                    "observed_eigenvalue": obs[j], "null_median": np.median(vals[:, j]),
                    "null_p95": np.quantile(vals[:, j], .95), "null_p99": np.quantile(vals[:, j], .99),
                    "empirical_upper_tail_p": (1 + np.sum(vals[:, j] >= obs[j])) / 10001,
                    "replicates": 10000,
                    "preserved": "axis marginal distributions" if name.startswith("INDEPENDENT") else "n,p,axis variance scale",
                    "broken": "cross-axis dependence" if name.startswith("INDEPENDENT") else "all empirical distribution and dependence structure",
                })
        p95a = np.quantile(nulls["INDEPENDENT_AXIS_PERMUTATION"], .95, axis=0)
        p95b = np.quantile(nulls["GAUSSIAN_REFERENCE"], .95, axis=0)
        thresholds[cohort] = {"a": p95a, "b": p95b, "observed": obs}
    return pd.DataFrame(rows), thresholds


def bootstrap_rank(resid: dict[str, np.ndarray], thresholds: dict[str, dict[str, np.ndarray]]) -> pd.DataFrame:
    rows: list[dict] = []
    for cohort in COHORTS:
        x = resid[cohort]
        n = len(x)
        rng = np.random.default_rng(910103 + (0 if cohort == "CGGA325" else 10000))
        _, u0 = eig(covariance(x))
        p95a, p95b = thresholds[cohort]["a"], thresholds[cohort]["b"]
        ranks = []
        cohort_rows = []
        for b in range(1, 5001):
            z = x[rng.integers(0, n, n)]
            values, vectors = eig(covariance(z))
            rank = effective_rank(values, p95a, p95b)
            ranks.append(rank)
            cohort_rows.append({
                "row_type": "REPLICATE", "cohort": cohort, "bootstrap_id": b, "effective_rank": rank,
                "rank_category": "RANK3_PLUS" if rank >= 3 else f"RANK{rank}",
                "lambda1": values[0], "lambda2": values[1], "lambda3": values[2],
                "lambda2_over_lambda3": values[1] / values[2],
                "rank2_projector_similarity": projection_similarity(u0, vectors, 2),
                "estimate": np.nan, "ci_low": np.nan, "ci_high": np.nan,
            })
        rows.extend(cohort_rows)
        ranks = np.asarray(ranks)
        for category, mask in [("RANK0", ranks == 0), ("RANK1", ranks == 1), ("RANK2", ranks == 2), ("RANK3_PLUS", ranks >= 3)]:
            rows.append({"row_type": "SUMMARY", "cohort": cohort, "bootstrap_id": np.nan, "effective_rank": np.nan,
                         "rank_category": category, "estimate": mask.mean(), "ci_low": np.nan, "ci_high": np.nan})
        for col in ["lambda1", "lambda2", "lambda3", "lambda2_over_lambda3", "rank2_projector_similarity"]:
            values = np.array([r[col] for r in cohort_rows])
            rows.append({"row_type": "SUMMARY", "cohort": cohort, "bootstrap_id": np.nan, "effective_rank": np.nan,
                         "rank_category": col.upper(), "estimate": np.median(values),
                         "ci_low": np.quantile(values, .025), "ci_high": np.quantile(values, .975)})
    return pd.DataFrame(rows)


def cv_rank_utility(resid: dict[str, np.ndarray]) -> pd.DataFrame:
    rows: list[dict] = []
    summary_rows: list[dict] = []
    for cohort in COHORTS:
        x = resid[cohort]
        n = len(x)
        rng = np.random.default_rng(910104 + (0 if cohort == "CGGA325" else 10000))
        patient_sse = {k: np.zeros(n) for k in (1, 2, 3)}
        patient_sst = np.zeros(n)
        patient_counts = np.zeros(n)
        for repeat in range(1, 101):
            all_pred = {k: np.empty_like(x) for k in (1, 2, 3)}
            baseline = np.empty_like(x)
            for test in folds(n, rng, 5):
                train = np.setdiff1d(np.arange(n), test, assume_unique=True)
                center, _, vectors = fit_basis(x[train])
                baseline[test] = center
                for k in (1, 2, 3):
                    all_pred[k][test] = reconstruct(x[test], center, vectors, k)
            sst = np.square(x - baseline)
            patient_sst += sst.sum(axis=1)
            patient_counts += 1
            for k in (1, 2, 3):
                se = np.square(x - all_pred[k])
                patient_sse[k] += se.sum(axis=1)
                q2, rmse = q2_rmse(x, all_pred[k], baseline)
                captured = 1 - np.trace(covariance(x - all_pred[k])) / np.trace(covariance(x - baseline))
                rows.append({"row_type": "REPEAT", "cohort": cohort, "repeat": repeat, "rank": k, "heldout_q2": q2,
                             "heldout_rmse": rmse, "heldout_covariance_captured": captured,
                             "delta_q2_rank2_minus_rank1": np.nan, "ci_low": np.nan, "ci_high": np.nan})
        delta_obs = (patient_sse[1].sum() - patient_sse[2].sum()) / patient_sst.sum()
        boot_rng = np.random.default_rng(910105 + (0 if cohort == "CGGA325" else 10000))
        bvals = np.empty(5000)
        for b in range(5000):
            ii = boot_rng.integers(0, n, n)
            bvals[b] = (patient_sse[1][ii].sum() - patient_sse[2][ii].sum()) / patient_sst[ii].sum()
        for k in (1, 2, 3):
            d = pd.DataFrame(rows)
            q = d[(d.row_type == "REPEAT") & (d.cohort == cohort) & (d["rank"] == k)]
            summary_rows.append({"row_type": "SUMMARY", "cohort": cohort, "repeat": np.nan, "rank": k,
                                 "heldout_q2": q.heldout_q2.mean(), "heldout_rmse": q.heldout_rmse.mean(),
                                 "heldout_covariance_captured": q.heldout_covariance_captured.mean(),
                                 "delta_q2_rank2_minus_rank1": delta_obs if k == 2 else np.nan,
                                 "ci_low": np.quantile(bvals, .025) if k == 2 else np.nan,
                                 "ci_high": np.quantile(bvals, .975) if k == 2 else np.nan})
    return pd.concat([pd.DataFrame(rows), pd.DataFrame(summary_rows)], ignore_index=True)


def cross_cohort_rank_transfer(resid: dict[str, np.ndarray]) -> pd.DataFrame:
    rows = []
    bases = {c: fit_basis(resid[c]) for c in COHORTS}
    for source, target in [("CGGA325", "CGGA693"), ("CGGA693", "CGGA325")]:
        center, values, vectors = bases[source]
        _, target_values, target_vectors = bases[target]
        y = resid[target]
        baseline = np.broadcast_to(center, y.shape)
        q2s = {}
        for k in (1, 2, 3):
            pred = reconstruct(y, center, vectors, k)
            q2, rmse = q2_rmse(y, pred, baseline)
            captured = 1 - np.trace(covariance(y - pred)) / np.trace(covariance(y - baseline))
            aa = angles_deg(vectors, target_vectors, k)
            q2s[k] = q2
            rows.append({
                "source_cohort": source, "target_cohort": target, "rank": k, "source_n": len(resid[source]), "target_n": len(y),
                "transferred_q2": q2, "transferred_rmse": rmse, "transferred_covariance_captured": captured,
                "principal_angles_deg": ";".join(f"{x:.9g}" for x in aa), "max_principal_angle_deg": aa.max(),
                "projection_similarity": projection_similarity(vectors, target_vectors, k),
                "grassmann_distance_rad": np.linalg.norm(np.radians(aa)),
                "delta_q2_rank2_minus_rank1": np.nan,
                "source_variance_fraction": values[:k].sum() / values.sum(),
                "target_native_variance_fraction": target_values[:k].sum() / target_values.sum(),
            })
        rows[-2]["delta_q2_rank2_minus_rank1"] = q2s[2] - q2s[1]
    return pd.DataFrame(rows)


def ppca_sensitivity(resid: dict[str, np.ndarray]) -> pd.DataFrame:
    rows = []
    for cohort, x in resid.items():
        n, p = x.shape
        values, _ = eig(covariance(x - x.mean(axis=0)))
        for q in range(0, p):
            if q == 0:
                sigma2 = values.mean()
                model_eigs = np.repeat(sigma2, p)
            else:
                sigma2 = values[q:].mean()
                model_eigs = np.r_[values[:q], np.repeat(sigma2, p - q)]
            logdet = np.log(model_eigs).sum()
            trace_term = np.sum(values / model_eigs)
            loglik = -0.5 * n * (p * np.log(2 * np.pi) + logdet + trace_term)
            npar = p * q - q * (q - 1) / 2 + 1
            bic = -2 * loglik + npar * np.log(n)
            rows.append({"cohort": cohort, "model": "PROBABILISTIC_PCA", "latent_rank": q, "noise_variance": sigma2,
                         "log_likelihood": loglik, "parameter_count": npar, "BIC": bic,
                         "delta_BIC_from_best": np.nan, "role": "GEOMETRIC_SENSITIVITY_ONLY"})
    out = pd.DataFrame(rows)
    out["delta_BIC_from_best"] = out.groupby("cohort").BIC.transform(lambda x: x - x.min())
    return out


def nearest_indices(query_s: np.ndarray, ref_s: np.ndarray, ref_ids: np.ndarray, exclude_self: bool = False) -> np.ndarray:
    out = np.empty(len(query_s), dtype=int)
    for i, value in enumerate(query_s):
        delta = np.abs(ref_s - value)
        if exclude_self:
            delta[i] = np.inf
        best = np.flatnonzero(delta == np.nanmin(delta))
        if len(best) > 1:
            out[i] = best[np.argmin(ref_ids[best])]
        else:
            out[i] = best[0]
    return out


def calibration_rows(master_by: dict[str, dict[str, np.ndarray]]) -> pd.DataFrame:
    rows = []
    for source, target in [("CGGA325", "CGGA693"), ("CGGA693", "CGGA325")]:
        ss, rs = master_by[source]["s"], master_by[source]["r"]
        st, rt = master_by[target]["s"], master_by[target]["r"]
        alpha, beta = slope_vector(ss, rs)
        pred = alpha + np.outer(st, beta)
        baseline = np.broadcast_to(alpha, rt.shape)
        q2, rmse = q2_rmse(rt, pred, baseline)
        y, p = rt.ravel(), pred.ravel()
        cal = np.linalg.lstsq(np.column_stack([np.ones(len(p)), p]), y, rcond=None)[0]
        rows.append({"row_type": "CROSS_COHORT_OVERALL", "source_cohort": source, "target_cohort": target, "axis": "ALL",
                     "n_source": len(ss), "n_target": len(st), "q2": q2, "rmse": rmse,
                     "intercept_only_rmse": math.sqrt(np.square(rt - baseline).mean()),
                     "calibration_intercept": cal[0], "calibration_slope": cal[1], "repeat": np.nan})
        for j, axis in enumerate(AXES):
            q2j, rmsej = q2_rmse(rt[:, [j]], pred[:, [j]], baseline[:, [j]])
            cj = np.linalg.lstsq(np.column_stack([np.ones(len(st)), pred[:, j]]), rt[:, j], rcond=None)[0]
            rows.append({"row_type": "CROSS_COHORT_AXIS", "source_cohort": source, "target_cohort": target, "axis": axis,
                         "n_source": len(ss), "n_target": len(st), "q2": q2j, "rmse": rmsej,
                         "intercept_only_rmse": math.sqrt(np.square(rt[:, j] - baseline[:, j]).mean()),
                         "calibration_intercept": cj[0], "calibration_slope": cj[1], "repeat": np.nan})
    for cohort in COHORTS:
        s, r = master_by[cohort]["s"], master_by[cohort]["r"]
        n = len(s)
        rng = np.random.default_rng(920105 + (0 if cohort == "CGGA325" else 10000))
        for repeat in range(1, 101):
            pred = np.empty_like(r)
            base = np.empty_like(r)
            for test in folds(n, rng, 5):
                train = np.setdiff1d(np.arange(n), test, assume_unique=True)
                alpha, beta = slope_vector(s[train], r[train])
                pred[test] = alpha + np.outer(s[test], beta)
                base[test] = alpha
            q2, rmse = q2_rmse(r, pred, base)
            rows.append({"row_type": "WITHIN_COHORT_CV", "source_cohort": cohort, "target_cohort": cohort, "axis": "ALL",
                         "n_source": n, "n_target": n, "q2": q2, "rmse": rmse,
                         "intercept_only_rmse": math.sqrt(np.square(r - base).mean()),
                         "calibration_intercept": np.nan, "calibration_slope": np.nan, "repeat": repeat})
    return pd.DataFrame(rows)


def r9b_analysis(master: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    data: dict[str, dict[str, np.ndarray]] = {}
    for cohort in COHORTS:
        d = master[master.cohort.eq(cohort)].sort_values("patient_id")
        data[cohort] = {"s": d.bulk_mes_z.to_numpy(float), "r": d[RAW_COLS].to_numpy(float), "ids": d.patient_id.astype(str).to_numpy()}
    beta = {}
    alpha = {}
    residual = {}
    sigma = {}
    for cohort in COHORTS:
        residual[cohort], alpha[cohort], beta[cohort] = residualize(data[cohort]["s"], data[cohort]["r"])
        sigma[cohort] = covariance(residual[cohort])

    obs_meaning = meaning_metrics(beta["CGGA325"], beta["CGGA693"])
    rng = np.random.default_rng(920101)
    beta_boot = {c: np.empty((5000, len(AXES))) for c in COHORTS}
    global_boot = {k: np.empty(5000) for k in obs_meaning}
    for b in range(5000):
        bb = {}
        for cohort in COHORTS:
            n = len(data[cohort]["s"])
            ii = rng.integers(0, n, n)
            bb[cohort] = slope_vector(data[cohort]["s"][ii], data[cohort]["r"][ii])[1]
            beta_boot[cohort][b] = bb[cohort]
        m = meaning_metrics(bb["CGGA325"], bb["CGGA693"])
        for k, v in m.items():
            global_boot[k][b] = v
    meaning_rows = []
    for cohort in COHORTS:
        for j, axis in enumerate(AXES):
            meaning_rows.append({"row_type": "AXIS_BETA", "cohort_or_contrast": cohort, "axis_or_metric": axis,
                                 "estimate": beta[cohort][j], "ci_low": np.quantile(beta_boot[cohort][:, j], .025),
                                 "ci_high": np.quantile(beta_boot[cohort][:, j], .975), "bootstrap_replicates": 5000,
                                 "empirical_p_vs_within_null": np.nan})
    for metric, value in obs_meaning.items():
        meaning_rows.append({"row_type": "BETWEEN_COHORT_GLOBAL", "cohort_or_contrast": "CGGA325_VS_CGGA693", "axis_or_metric": metric,
                             "estimate": value, "ci_low": np.quantile(global_boot[metric], .025),
                             "ci_high": np.quantile(global_boot[metric], .975), "bootstrap_replicates": 5000,
                             "empirical_p_vs_within_null": np.nan})

    split_rows = []
    rng = np.random.default_rng(920102)
    pooled_null = {k: [] for k in obs_meaning}
    for cohort in COHORTS:
        s, r = data[cohort]["s"], data[cohort]["r"]
        n = len(s)
        for b in range(1, 5001):
            order = rng.permutation(n)
            a, z = np.array_split(order, 2)
            ba = slope_vector(s[a], r[a])[1]
            bz = slope_vector(s[z], r[z])[1]
            m = meaning_metrics(ba, bz)
            split_rows.append({"cohort": cohort, "split_id": b, "n_half_a": len(a), "n_half_b": len(z), **m})
            for k, v in m.items():
                pooled_null[k].append(v)
    for row in meaning_rows:
        if row["row_type"] != "BETWEEN_COHORT_GLOBAL":
            continue
        metric, value = row["axis_or_metric"], row["estimate"]
        null = np.asarray(pooled_null[metric])
        if "similarity" in metric or "congruence" in metric:
            p = (1 + np.sum(null <= value)) / (len(null) + 1)
        else:
            p = (1 + np.sum(null >= value)) / (len(null) + 1)
        row["empirical_p_vs_within_null"] = p

    obs_geom = geometry_metrics(sigma["CGGA325"], sigma["CGGA693"])
    geom_boot = {k: np.empty(5000) for k in obs_geom}
    rng = np.random.default_rng(920103)
    for b in range(5000):
        covs = {}
        for cohort in COHORTS:
            x = residual[cohort]
            covs[cohort] = covariance(x[rng.integers(0, len(x), len(x))])
        m = geometry_metrics(covs["CGGA325"], covs["CGGA693"])
        for k, v in m.items():
            geom_boot[k][b] = v
    pool = np.vstack([residual[c] for c in COHORTS])
    n325 = len(residual["CGGA325"])
    geom_null = {k: np.empty(5000) for k in obs_geom}
    for b in range(5000):
        order = rng.permutation(len(pool))
        m = geometry_metrics(covariance(pool[order[:n325]]), covariance(pool[order[n325:]]))
        for k, v in m.items():
            geom_null[k][b] = v
    geom_rows = []
    for metric, value in obs_geom.items():
        similarity_metric = metric in {"rv_coefficient", "matrix_correlation_upper", "rank2_projection_similarity"}
        p = (1 + np.sum(geom_null[metric] <= value)) / 5001 if similarity_metric else (1 + np.sum(geom_null[metric] >= value)) / 5001
        geom_rows.append({"metric": metric, "estimate": value, "bootstrap_ci_low": np.quantile(geom_boot[metric], .025),
                          "bootstrap_ci_high": np.quantile(geom_boot[metric], .975), "permutation_null_median": np.median(geom_null[metric]),
                          "permutation_null_p05": np.quantile(geom_null[metric], .05), "permutation_null_p95": np.quantile(geom_null[metric], .95),
                          "empirical_noninvariance_p": p, "bootstrap_replicates": 5000, "permutation_replicates": 5000})

    nearest_rows = []
    nearest_summary = []
    for source, target in [("CGGA325", "CGGA693"), ("CGGA693", "CGGA325")]:
        q = data[source]
        t = data[target]
        cross_idx = nearest_indices(q["s"], t["s"], t["ids"], False)
        same_idx = nearest_indices(q["s"], q["s"], q["ids"], True)
        cross_dist = np.linalg.norm(q["r"] - t["r"][cross_idx], axis=1)
        same_dist = np.linalg.norm(q["r"] - q["r"][same_idx], axis=1)
        delta = cross_dist - same_dist
        for i in range(len(q["s"])):
            nearest_rows.append({"row_type": "PATIENT", "direction": f"{source}_TO_{target}", "source_patient_id": q["ids"][i],
                                 "source_mes_z": q["s"][i], "cross_match_patient_id": t["ids"][cross_idx[i]],
                                 "cross_delta_mes": abs(q["s"][i] - t["s"][cross_idx[i]]), "cross_realization_distance": cross_dist[i],
                                 "same_match_patient_id": q["ids"][same_idx[i]], "same_delta_mes": abs(q["s"][i] - q["s"][same_idx[i]]),
                                 "same_realization_distance": same_dist[i], "cross_minus_same_distance": delta[i],
                                 "estimate": np.nan, "ci_low": np.nan, "ci_high": np.nan, "empirical_p": np.nan})
        rngn = np.random.default_rng(920104 + (0 if source == "CGGA325" else 10000))
        bv = np.empty(5000)
        for b in range(5000):
            ii = rngn.integers(0, len(delta), len(delta))
            bv[b] = np.median(delta[ii])
        nearest_summary.append({"row_type": "SUMMARY", "direction": f"{source}_TO_{target}", "source_patient_id": "PATIENT_MEDIAN",
                                "source_mes_z": np.nan, "cross_match_patient_id": "NA", "cross_delta_mes": np.median(np.abs(q["s"] - t["s"][cross_idx])),
                                "cross_realization_distance": np.median(cross_dist), "same_match_patient_id": "NA",
                                "same_delta_mes": np.median(np.abs(q["s"] - q["s"][same_idx])), "same_realization_distance": np.median(same_dist),
                                "cross_minus_same_distance": np.median(delta), "estimate": np.median(delta),
                                "ci_low": np.quantile(bv, .025), "ci_high": np.quantile(bv, .975),
                                "empirical_p": 2 * min((1 + np.sum(bv <= 0)) / 5001, (1 + np.sum(bv >= 0)) / 5001)})

    calibration = calibration_rows(data)
    diagnostics = {"beta": {k: beta[k].tolist() for k in COHORTS}, "meaning": obs_meaning, "geometry": obs_geom,
                   "nearest_summary": nearest_summary}
    return (pd.DataFrame(meaning_rows), pd.DataFrame(split_rows), pd.DataFrame(geom_rows),
            pd.concat([pd.DataFrame(nearest_rows), pd.DataFrame(nearest_summary)], ignore_index=True), calibration, diagnostics)


def main() -> None:
    AOUT.mkdir(parents=True, exist_ok=True)
    BOUT.mkdir(parents=True, exist_ok=True)
    LOG.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(PATIENT_REALIZATION / "01_reference/PATIENT_REALIZATION_MASTER.tsv", sep="\t", dtype={"patient_id": str})
    fp = pd.read_csv(PATIENT_REALIZATION / "02_fingerprint/MES_CONDITIONED_FINGERPRINT.tsv", sep="\t", dtype={"patient_id": str})
    if master.duplicated(["cohort", "patient_id"]).any() or fp.duplicated(["cohort", "patient_id"]).any():
        raise RuntimeError("Patient key uniqueness check did not pass")
    if master[RAW_COLS + ["bulk_mes_z"]].isna().any().any() or fp[RESID_COLS].isna().any().any():
        raise RuntimeError("Missing primary R9A/R9B value")
    resid = {c: fp.loc[fp.cohort.eq(c)].sort_values("patient_id")[RESID_COLS].to_numpy(float) for c in COHORTS}

    pa, thresholds = pa_analysis(resid)
    write(pa, AOUT / "R9A_PARALLEL_ANALYSIS.tsv")
    boot = bootstrap_rank(resid, thresholds)
    write(boot, AOUT / "R9A_BOOTSTRAP_RANK_DISTRIBUTION.tsv")
    cv = cv_rank_utility(resid)
    write(cv, AOUT / "R9A_CV_RANK_UTILITY.tsv")
    transfer = cross_cohort_rank_transfer(resid)
    write(transfer, AOUT / "R9A_CROSS_COHORT_RANK_TRANSFER.tsv")
    perturb = pd.DataFrame([{
        "status": "NOT_EVALUABLE_GENE_PERTURBATION", "scope": "COMPLETE_FIVE_AXIS_PERTURBATION",
        "reason": "MYELOID and VASCULAR_STROMAL are frozen cell-fraction axes without gene-set memberships; perturbing only score-derived axes would not satisfy the frozen complete-five-axis test.",
        "conditions_requested": "retain_90_percent;retain_80_percent;leave_small_block_out", "replicates_run": 0,
        "analysis_decision": "NO_PARTIAL_SURROGATE_ANALYSIS",
    }])
    write(perturb, AOUT / "R9A_GENE_PERTURBATION_STABILITY.tsv")
    model = ppca_sensitivity(resid)
    write(model, AOUT / "R9A_MODEL_SENSITIVITY.tsv")

    rank_obs = {c: effective_rank(thresholds[c]["observed"], thresholds[c]["a"], thresholds[c]["b"]) for c in COHORTS}
    rank2_prob = {c: float(boot[(boot.row_type == "SUMMARY") & (boot.cohort == c) & (boot.rank_category == "RANK2")].estimate.iloc[0]) for c in COHORTS}
    cvsum = cv[(cv.row_type == "SUMMARY") & (cv["rank"] == 2)].set_index("cohort")
    t2 = transfer[transfer["rank"] == 2].set_index("source_cohort")
    strong = all(rank_obs[c] == 2 for c in COHORTS) and all(rank2_prob[c] >= .80 for c in COHORTS) and bool((cvsum.ci_low > 0).all()) and bool((t2.delta_q2_rank2_minus_rank1 > 0).all())
    adjudication = "LOW_DIMENSIONAL_REALIZATION_SUBSPACE_STRONGLY_SUPPORTED" if strong else "LOW_DIMENSIONAL_STRUCTURE_SUPPORTED_EFFECTIVE_RANK_UNCERTAIN"
    (AOUT / "R9A_STRUCTURAL_VALIDITY_ADJUDICATION.md").write_text(
        "# R9A structural-validity adjudication\n\n"
        f"Status: `{adjudication}`\n\n"
        f"- Dual-null effective rank: CGGA325={rank_obs['CGGA325']}; CGGA693={rank_obs['CGGA693']}.\n"
        f"- Bootstrap rank-2 probability: CGGA325={rank2_prob['CGGA325']:.3f}; CGGA693={rank2_prob['CGGA693']:.3f}.\n"
        f"- Held-out delta-Q2 rank2-rank1 (95% patient-bootstrap CI): CGGA325={cvsum.loc['CGGA325','delta_q2_rank2_minus_rank1']:.4f} [{cvsum.loc['CGGA325','ci_low']:.4f}, {cvsum.loc['CGGA325','ci_high']:.4f}]; CGGA693={cvsum.loc['CGGA693','delta_q2_rank2_minus_rank1']:.4f} [{cvsum.loc['CGGA693','ci_low']:.4f}, {cvsum.loc['CGGA693','ci_high']:.4f}].\n"
        f"- Bidirectional transferred delta-Q2 rank2-rank1: CGGA325→CGGA693={t2.loc['CGGA325','delta_q2_rank2_minus_rank1']:.4f}; CGGA693→CGGA325={t2.loc['CGGA693','delta_q2_rank2_minus_rank1']:.4f}.\n"
        "- Complete five-axis gene-set perturbation is not evaluable because two frozen axes are cell fractions without gene memberships; no partial surrogate was substituted.\n"
        "- R8 external structural evidence remains `PENDING_R8_EXTERNAL_INTEGRATION` unless its final QA gate closes.\n\n"
        "Allowed wording: canonical MES behaves as a scalar projection of a reproducible low-dimensional biological realization space. Prohibited wording: MES itself is rank-2 or defines two subtypes.\n",
        encoding="utf-8",
    )

    meaning, split, geom, nearest, cal, diag = r9b_analysis(master)
    write(meaning, BOUT / "R9B_MEANING_VECTOR_INVARIANCE.tsv")
    write(split, BOUT / "R9B_WITHIN_COHORT_REPRODUCIBILITY_NULL.tsv")
    write(geom, BOUT / "R9B_RESIDUAL_GEOMETRY_INVARIANCE.tsv")
    write(nearest, BOUT / "R9B_CONDITIONAL_REALIZATION_TRANSPORT.tsv")
    write(cal, BOUT / "R9B_CROSS_COHORT_BIOLOGICAL_CALIBRATION.tsv")

    mdist = meaning[(meaning.row_type == "BETWEEN_COHORT_GLOBAL") & (meaning.axis_or_metric == "standardized_euclidean_distance")].iloc[0]
    cosrow = meaning[(meaning.row_type == "BETWEEN_COHORT_GLOBAL") & (meaning.axis_or_metric == "cosine_similarity")].iloc[0]
    beta_direction = cosrow.estimate >= .90
    beta_within = mdist.empirical_p_vs_within_null >= .05
    geom_p = geom.set_index("metric").loc["trace_normalized_frobenius_distance", "empirical_noninvariance_p"]
    residual_invariant = geom_p >= .05
    ns = nearest[nearest.row_type == "SUMMARY"]
    cross_greater = bool((ns.ci_low > 0).all())
    if beta_direction and beta_within and residual_invariant:
        bstatus = "BIOLOGICAL_MEANING_INVARIANT_WITHIN_REPRODUCIBILITY_LIMIT"
    elif beta_direction and not residual_invariant:
        bstatus = "DIRECTIONAL_MEANING_INVARIANT_RESIDUAL_GEOMETRY_NONINVARIANT"
    elif beta_direction:
        bstatus = "PARTIAL_MEASUREMENT_INVARIANCE"
    else:
        bstatus = "BIOLOGICAL_MEANING_NONINVARIANT"
    (BOUT / "R9B_MEASUREMENT_INVARIANCE_ADJUDICATION.md").write_text(
        "# R9B measurement-invariance adjudication\n\n"
        f"Status: `{bstatus}`\n\n"
        f"- Meaning-vector cosine similarity={cosrow.estimate:.4f} (95% bootstrap CI {cosrow.ci_low:.4f} to {cosrow.ci_high:.4f}).\n"
        f"- Between-cohort standardized Euclidean beta distance={mdist.estimate:.4f}; empirical comparison with pooled within-cohort split-half distances P={mdist.empirical_p_vs_within_null:.4g}.\n"
        f"- Residual covariance trace-normalized Frobenius distance={diag['geometry']['trace_normalized_frobenius_distance']:.4f}; permutation non-invariance P={geom_p:.4g}.\n"
        f"- Cross-cohort nearest-MES realization distance exceeded the same-cohort reference in both directions: {'YES' if cross_greater else 'NO'}.\n"
        "- A transported directional beta vector does not imply biological identity: conditional realization freedom and calibration are adjudicated separately.\n",
        encoding="utf-8",
    )

    env = {
        "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
        "platform": platform.platform(), "R9A_adjudication": adjudication, "R9B_adjudication": bstatus,
        "visualization": "NOT_REQUESTED_BY_USER",
    }
    (LOG / "10_r9a_r9b_environment.json").write_text(json.dumps(env, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"R9A": adjudication, "R9B": bstatus, "rank": rank_obs, "rank2_probability": rank2_prob}, ensure_ascii=False))


if __name__ == "__main__":
    main()
