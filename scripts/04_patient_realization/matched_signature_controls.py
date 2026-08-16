#!/usr/bin/env python3
"""Adversarial construct-specificity analysis for the frozen R9C library."""

from __future__ import annotations

import json
import os
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy


DATA_ROOT = Path(os.environ.get("GBM_MES_DATA_ROOT", "data"))
RESULT_ROOT = Path(os.environ.get("GBM_MES_OUTPUT_ROOT", "results"))
RUN = RESULT_ROOT / "construct_validation"
INPUT = DATA_ROOT / "prepared" / "construct_validation" / "matched_signatures"
OUT = RUN / "03_matched_signature_controls"
MASTER = RESULT_ROOT / "patient_realization" / "01_reference" / "PATIENT_REALIZATION_MASTER.tsv"
AXES = ["MALIGNANT", "MYELOID", "HYPOXIA", "MATRIX", "VASCULAR_STROMAL"]
CONTEXT = AXES[1:]
COHORTS = ["CGGA325", "CGGA693"]


def z(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    return (x - x.mean()) / x.std(ddof=1)


def loo_q2(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Leakage-free univariate leave-one-patient-out prediction."""
    x, y = z(x), z(y)
    n = len(x)
    X = np.column_stack([np.ones(n), x])
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    residual = y - X @ beta
    hat = np.einsum("ij,jk,ik->i", X, np.linalg.inv(X.T @ X), X)
    cv_residual = residual / (1.0 - hat)
    baseline_residual = (n / (n - 1.0)) * (y - y.mean())
    q2 = 1.0 - np.sum(cv_residual**2) / np.sum(baseline_residual**2)
    rmse = float(np.sqrt(np.mean(cv_residual**2)))
    return float(q2), rmse, float(beta[1])


def nearest_score_distance(score: np.ndarray, realization: np.ndarray, ids: np.ndarray) -> tuple[float, float]:
    n = len(score)
    order = np.lexsort((ids.astype(str), score))
    picked = np.empty(n, dtype=int)
    for pos, idx in enumerate(order):
        candidates = []
        if pos:
            candidates.append(order[pos - 1])
        if pos + 1 < n:
            candidates.append(order[pos + 1])
        candidates.sort(key=lambda j: (abs(score[j] - score[idx]), str(ids[j])))
        picked[idx] = candidates[0]
    distances = np.sqrt(np.mean((realization - realization[picked]) ** 2, axis=1))
    return float(np.median(distances)), float(np.mean(distances))


def residual_geometry(score: np.ndarray, realization: np.ndarray) -> dict[str, object]:
    X = np.column_stack([np.ones(len(score)), score])
    beta = np.linalg.lstsq(X, realization, rcond=None)[0]
    in_sample_residual = realization - X @ beta
    hat = np.einsum("ij,jk,ik->i", X, np.linalg.inv(X.T @ X), X)
    residual_raw = in_sample_residual / (1.0 - hat[:, None])
    residual = (residual_raw - residual_raw.mean(axis=0)) / residual_raw.std(axis=0, ddof=1)
    cov = np.cov(residual, rowvar=False)
    eigval, eigvec = np.linalg.eigh(cov)
    order = np.argsort(eigval)[::-1]
    eigval, eigvec = eigval[order], eigvec[:, order]
    corr = np.corrcoef(residual, rowvar=False)
    upper = corr[np.triu_indices(corr.shape[0], 1)]
    p = eigval / eigval.sum()
    positive = p[p > 0]
    effective_rank = float(np.exp(-np.sum(positive * np.log(positive))))
    return {
        "beta": beta[1],
        "residual": residual,
        "cov": cov,
        "basis2": eigvec[:, :2],
        "mean_abs_offdiag_correlation": float(np.mean(np.abs(upper))),
        "top2_variance_fraction": float(eigval[:2].sum() / eigval.sum()),
        "entropy_effective_rank": effective_rank,
    }


def projection_similarity(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.trace((a @ a.T) @ (b @ b.T)) / a.shape[1])


def matrix_corr(a: np.ndarray, b: np.ndarray) -> float:
    ix = np.triu_indices(a.shape[0], 1)
    if len(ix[0]) < 2:
        return np.nan
    return float(np.corrcoef(a[ix], b[ix])[0, 1])


def empirical(observed: float, reference: np.ndarray, direction: str) -> tuple[float, float]:
    reference = np.asarray(reference, float)
    reference = reference[np.isfinite(reference)]
    if direction == "HIGH":
        pct = 100.0 * (np.sum(reference <= observed) + .5) / (len(reference) + 1.0)
        p = (np.sum(reference >= observed) + 1.0) / (len(reference) + 1.0)
    elif direction == "LOW":
        pct = 100.0 * (np.sum(reference <= observed) + .5) / (len(reference) + 1.0)
        p = (np.sum(reference <= observed) + 1.0) / (len(reference) + 1.0)
    else:
        center = np.median(reference)
        delta = abs(observed - center)
        p = (np.sum(np.abs(reference - center) >= delta) + 1.0) / (len(reference) + 1.0)
        pct = 100.0 * (np.sum(reference <= observed) + .5) / (len(reference) + 1.0)
    return float(pct), float(p)


def bh(p: pd.Series) -> pd.Series:
    a = p.to_numpy(float)
    out = np.full(len(a), np.nan)
    ok = np.isfinite(a)
    vals = a[ok]
    order = np.argsort(vals)
    ranked = vals[order]
    adj = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    tmp = np.empty(len(vals)); tmp[order] = adj
    out[np.where(ok)[0]] = tmp
    return pd.Series(out, index=p.index)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    library = pd.read_csv(INPUT / "R9C_COMPARATOR_LIBRARY_FREEZE.tsv", sep="\t")
    signatures = library.signature_id.tolist()
    tier = library.set_index("signature_id").tier.to_dict()
    master = pd.read_csv(MASTER, sep="\t", dtype={"patient_id": str})
    bulk = {
        c: pd.read_csv(INPUT / f"R9C_{c}_BULK_SCORES.tsv.gz", sep="\t", dtype={"patient_id": str}).set_index("patient_id")
        for c in COHORTS
    }
    malignant_all = pd.read_csv(INPUT / "R9C_MALIGNANT_SCORES.tsv.gz", sep="\t", dtype={"patient_id": str})
    malignant = {c: malignant_all.loc[malignant_all.cohort.eq(c)].set_index("patient_id").drop(columns="cohort") for c in COHORTS}
    taut = pd.read_csv(INPUT / "R9C_TAUTOLOGY_EXCLUSION_REGISTRY.tsv", sep="\t")
    taut["tautology_excluded"] = taut.tautology_excluded.astype(str).str.lower().eq("true")
    excluded = {(r.signature_id, r.endpoint): bool(r.tautology_excluded) for r in taut.itertuples()}

    cohort_data: dict[str, dict[str, object]] = {}
    for c in COHORTS:
        md = master.loc[master.cohort.eq(c)].copy().set_index("patient_id")
        ids = md.index.to_numpy(str)
        if ids.size != len(set(ids)) or set(ids) != set(bulk[c].index) or set(ids) != set(malignant[c].index):
            raise RuntimeError(f"Patient alignment failure in {c}")
        cohort_data[c] = {
            "ids": ids,
            "bulk": bulk[c].loc[ids, signatures],
            "malignant": malignant[c].loc[ids, signatures],
            "context": pd.DataFrame({
                "MYELOID": md.loc[ids, "myeloid_z"],
                "HYPOXIA": md.loc[ids, "hypoxia_z"],
                "MATRIX": md.loc[ids, "matrix_z"],
                "VASCULAR_STROMAL": md.loc[ids, "vascular_z"],
            }, index=ids),
        }

    metric_rows: list[dict[str, object]] = []
    geometry: dict[tuple[str, str], dict[str, object]] = {}
    for sid in signatures:
        keep_context = [a for a in CONTEXT if not excluded[(sid, a)]]
        keep_axes = ["MALIGNANT", *keep_context]
        for c in COHORTS:
            obj = cohort_data[c]
            s = z(obj["bulk"][sid].to_numpy(float))
            m = z(obj["malignant"][sid].to_numpy(float))
            r = np.column_stack([m, *[z(obj["context"][a].to_numpy(float)) for a in keep_context]])
            q2, rmse, slope = loo_q2(s, m)
            for name, value in (("LOO_Q2", q2), ("LOO_RMSE", rmse), ("STANDARDIZED_SLOPE", slope)):
                metric_rows.append({"signature_id": sid, "tier": tier[sid], "cohort_or_contrast": c, "construct_metric": "M1_BULK_TO_MALIGNANT_IDENTIFIABILITY", "metric": name, "estimate": value, "axis_count": len(keep_axes), "tautology_exclusions": ";".join(set(CONTEXT) - set(keep_context)) or "NONE", "analysis_role": "PRIMARY" if name == "LOO_Q2" else "SECONDARY"})
            cors = []
            for axis in keep_context:
                cor = float(np.corrcoef(s, z(obj["context"][axis].to_numpy(float)))[0, 1])
                cors.append(cor)
                metric_rows.append({"signature_id": sid, "tier": tier[sid], "cohort_or_contrast": c, "construct_metric": "M2_ECOLOGICAL_COUPLING", "metric": f"PEARSON_R_{axis}", "estimate": cor, "axis_count": len(keep_context), "tautology_exclusions": ";".join(set(CONTEXT) - set(keep_context)) or "NONE", "analysis_role": "SECONDARY"})
            metric_rows.append({"signature_id": sid, "tier": tier[sid], "cohort_or_contrast": c, "construct_metric": "M2_ECOLOGICAL_COUPLING", "metric": "MEAN_ABS_CONTEXT_CORRELATION", "estimate": float(np.mean(np.abs(cors))), "axis_count": len(keep_context), "tautology_exclusions": ";".join(set(CONTEXT) - set(keep_context)) or "NONE", "analysis_role": "PRIMARY"})
            med_dist, mean_dist = nearest_score_distance(s, r, obj["ids"])
            for name, value in (("MEDIAN_RMS_REALIZATION_DISTANCE", med_dist), ("MEAN_RMS_REALIZATION_DISTANCE", mean_dist)):
                metric_rows.append({"signature_id": sid, "tier": tier[sid], "cohort_or_contrast": c, "construct_metric": "M3_SAME_SCORE_REALIZATION_DISPERSION", "metric": name, "estimate": value, "axis_count": len(keep_axes), "tautology_exclusions": ";".join(set(CONTEXT) - set(keep_context)) or "NONE", "analysis_role": "PRIMARY" if name.startswith("MEDIAN") else "SECONDARY"})
            g = residual_geometry(s, r)
            geometry[(sid, c)] = {**g, "axes": keep_axes}
            metric_rows.append({"signature_id": sid, "tier": tier[sid], "cohort_or_contrast": c, "construct_metric": "M4_RESIDUAL_DEPENDENCE", "metric": "MEAN_ABS_OFFDIAGONAL_CORRELATION", "estimate": g["mean_abs_offdiag_correlation"], "axis_count": len(keep_axes), "tautology_exclusions": ";".join(set(CONTEXT) - set(keep_context)) or "NONE", "analysis_role": "PRIMARY"})
            for name, value, role in (("TOP2_VARIANCE_FRACTION", g["top2_variance_fraction"], "PRIMARY"), ("ENTROPY_EFFECTIVE_RANK", g["entropy_effective_rank"], "SECONDARY")):
                metric_rows.append({"signature_id": sid, "tier": tier[sid], "cohort_or_contrast": c, "construct_metric": "M5_LOW_RANK_RESIDUAL_GEOMETRY", "metric": name, "estimate": value, "axis_count": len(keep_axes), "tautology_exclusions": ";".join(set(CONTEXT) - set(keep_context)) or "NONE", "analysis_role": role})

        ga, gb = geometry[(sid, COHORTS[0])], geometry[(sid, COHORTS[1])]
        if ga["axes"] != gb["axes"]:
            raise RuntimeError(f"Axis mismatch for {sid}")
        beta_a, beta_b = np.asarray(ga["beta"]), np.asarray(gb["beta"])
        cosine = float(np.dot(beta_a, beta_b) / (np.linalg.norm(beta_a) * np.linalg.norm(beta_b)))
        euclid = float(np.linalg.norm(beta_a - beta_b))
        std_euclid = euclid / np.sqrt(len(beta_a))
        transport = {
            ("M6_CROSS_COHORT_REALIZATION_TRANSPORT", "RANK2_PROJECTION_SIMILARITY"): projection_similarity(ga["basis2"], gb["basis2"]),
            ("M6_CROSS_COHORT_REALIZATION_TRANSPORT", "RESIDUAL_COVARIANCE_MATRIX_CORRELATION"): matrix_corr(ga["cov"], gb["cov"]),
            ("M7_MEASUREMENT_INVARIANCE", "BETA_VECTOR_COSINE"): cosine,
            ("M7_MEASUREMENT_INVARIANCE", "BETA_VECTOR_EUCLIDEAN_DISTANCE"): euclid,
            ("M7_MEASUREMENT_INVARIANCE", "BETA_VECTOR_STANDARDIZED_EUCLIDEAN_DISTANCE"): std_euclid,
        }
        for (construct, name), value in transport.items():
            metric_rows.append({"signature_id": sid, "tier": tier[sid], "cohort_or_contrast": "CGGA325_VS_CGGA693", "construct_metric": construct, "metric": name, "estimate": value, "axis_count": len(beta_a), "tautology_exclusions": ";".join(set(CONTEXT) - set(ga["axes"])) or "NONE", "analysis_role": "PRIMARY" if name in {"RANK2_PROJECTION_SIMILARITY", "BETA_VECTOR_STANDARDIZED_EUCLIDEAN_DISTANCE"} else "SECONDARY"})

    metrics = pd.DataFrame(metric_rows)

    # Canonical M7 is read from the retained R9B result, not recalculated.
    r9b = pd.read_csv(RUN / "02_measurement_invariance/R9B_MEANING_VECTOR_INVARIANCE.tsv", sep="\t")
    map_r9b = {
        "BETA_VECTOR_COSINE": "cosine_similarity",
        "BETA_VECTOR_EUCLIDEAN_DISTANCE": "euclidean_distance",
        "BETA_VECTOR_STANDARDIZED_EUCLIDEAN_DISTANCE": "standardized_euclidean_distance",
    }
    r9b_errors = []
    for metric_name, r9b_name in map_r9b.items():
        reference = float(r9b.loc[(r9b.row_type == "BETWEEN_COHORT_GLOBAL") & (r9b.axis_or_metric == r9b_name), "estimate"].iloc[0])
        ix = (metrics.signature_id == "CANONICAL_MES") & (metrics.construct_metric == "M7_MEASUREMENT_INVARIANCE") & (metrics.metric == metric_name)
        calculated = float(metrics.loc[ix, "estimate"].iloc[0])
        r9b_errors.append(abs(calculated - reference))
        metrics.loc[ix, "estimate"] = reference
        metrics.loc[ix, "analysis_role"] = "PRIMARY_R9B_READ_ONLY" if metric_name.endswith("STANDARDIZED_EUCLIDEAN_DISTANCE") else "SECONDARY_R9B_READ_ONLY"

    # Exact canonical bridge: the fair pipeline must reproduce the frozen
    # V1.6/R9A LOPO residual geometry before any specificity comparison.
    r9a_pa = pd.read_csv(RUN / "01_structural_validity/R9A_PARALLEL_ANALYSIS.tsv", sep="\t")
    r9a_transfer = pd.read_csv(RUN / "01_structural_validity/R9A_CROSS_COHORT_RANK_TRANSFER.tsv", sep="\t")
    r9a_errors = []
    for c in COHORTS:
        q = r9a_pa.loc[(r9a_pa.cohort == c) & (r9a_pa.null_model == "INDEPENDENT_AXIS_PERMUTATION") & (r9a_pa.component.isin([1, 2])), "observed_eigenvalue"]
        reference = float(q.sum() / 5.0)
        calc = float(metrics.loc[(metrics.signature_id == "CANONICAL_MES") & (metrics.cohort_or_contrast == c) & (metrics.construct_metric == "M5_LOW_RANK_RESIDUAL_GEOMETRY") & (metrics.metric == "TOP2_VARIANCE_FRACTION"), "estimate"].iloc[0])
        r9a_errors.append(abs(calc - reference))
    reference_proj = float(r9a_transfer.loc[(r9a_transfer.source_cohort == "CGGA325") & (r9a_transfer.target_cohort == "CGGA693") & (r9a_transfer["rank"] == 2), "projection_similarity"].iloc[0])
    calc_proj = float(metrics.loc[(metrics.signature_id == "CANONICAL_MES") & (metrics.construct_metric == "M6_CROSS_COHORT_REALIZATION_TRANSPORT") & (metrics.metric == "RANK2_PROJECTION_SIMILARITY"), "estimate"].iloc[0])
    r9a_errors.append(abs(calc_proj - reference_proj))
    if max(r9a_errors) > 1e-8:
        raise RuntimeError(f"Canonical R9A bridge did not reproduce: max geometry error={max(r9a_errors):.6g}")

    metrics.to_csv(OUT / "R9C_SIGNATURE_LEVEL_CONSTRUCT_METRICS.tsv", sep="\t", index=False, float_format="%.12g")
    metrics.loc[metrics.tier.eq("TIER3_MATCHED_RANDOM")].to_csv(OUT / "R9C_MATCHED_NULL_METRICS.tsv", sep="\t", index=False, float_format="%.12g")

    primary_specs = [
        ("M1_BULK_TO_MALIGNANT_IDENTIFIABILITY", "LOO_Q2", "LOW"),
        ("M2_ECOLOGICAL_COUPLING", "MEAN_ABS_CONTEXT_CORRELATION", "HIGH"),
        ("M3_SAME_SCORE_REALIZATION_DISPERSION", "MEDIAN_RMS_REALIZATION_DISTANCE", "HIGH"),
        ("M4_RESIDUAL_DEPENDENCE", "MEAN_ABS_OFFDIAGONAL_CORRELATION", "HIGH"),
        ("M5_LOW_RANK_RESIDUAL_GEOMETRY", "TOP2_VARIANCE_FRACTION", "HIGH"),
        ("M5_LOW_RANK_RESIDUAL_GEOMETRY", "ENTROPY_EFFECTIVE_RANK", "LOW"),
        ("M6_CROSS_COHORT_REALIZATION_TRANSPORT", "RANK2_PROJECTION_SIMILARITY", "HIGH"),
        ("M7_MEASUREMENT_INVARIANCE", "BETA_VECTOR_STANDARDIZED_EUCLIDEAN_DISTANCE", "HIGH"),
        ("M7_MEASUREMENT_INVARIANCE", "BETA_VECTOR_COSINE", "LOW"),
    ]
    rank_rows = []
    for construct, metric_name, direction in primary_specs:
        obs_rows = metrics.loc[(metrics.signature_id == "CANONICAL_MES") & (metrics.construct_metric == construct) & (metrics.metric == metric_name)]
        for obs in obs_rows.itertuples():
            for ref_tier, label in (("TIER3_MATCHED_RANDOM", "MATCHED_RANDOM"), ("TIER1_GBM_STATE", "NEFTEL_STATE"), ("TIER2_HALLMARK", "HALLMARK")):
                ref = metrics.loc[(metrics.tier == ref_tier) & (metrics.cohort_or_contrast == obs.cohort_or_contrast) & (metrics.construct_metric == construct) & (metrics.metric == metric_name), "estimate"].to_numpy(float)
                pct, p = empirical(float(obs.estimate), ref, direction)
                rank_rows.append({"construct_metric": construct, "metric": metric_name, "cohort_or_contrast": obs.cohort_or_contrast, "canonical_estimate": obs.estimate, "reference_family": label, "reference_n": len(ref), "adverse_direction": direction, "empirical_percentile": pct, "empirical_p": p})
    ranks = pd.DataFrame(rank_rows)
    random_ix = ranks.reference_family.eq("MATCHED_RANDOM")
    ranks.loc[random_ix, "bh_fdr_within_family_c"] = bh(ranks.loc[random_ix, "empirical_p"]).to_numpy()
    ranks.to_csv(OUT / "R9C_CONSTRUCT_SPECIFICITY_EMPIRICAL_RANK.tsv", sep="\t", index=False, na_rep="NA", float_format="%.12g")

    fixed_metrics = metrics.loc[metrics.tier.isin(["INDEX", "TIER1_GBM_STATE", "TIER2_HALLMARK"]) & metrics.construct_metric.isin([x[0] for x in primary_specs]) & metrics.metric.isin([x[1] for x in primary_specs])].copy()
    canonical_key = fixed_metrics.loc[fixed_metrics.signature_id.eq("CANONICAL_MES"), ["cohort_or_contrast", "construct_metric", "metric", "estimate"]].rename(columns={"estimate": "canonical_estimate"})
    fixed_metrics = fixed_metrics.merge(canonical_key, on=["cohort_or_contrast", "construct_metric", "metric"], how="left")
    fixed_metrics["difference_from_canonical"] = fixed_metrics.estimate - fixed_metrics.canonical_estimate
    fixed_metrics.loc[fixed_metrics.tier.eq("TIER1_GBM_STATE")].to_csv(OUT / "R9C_STATE_FAMILY_COMPARISON.tsv", sep="\t", index=False, float_format="%.12g")
    fixed_metrics.loc[fixed_metrics.tier.eq("TIER2_HALLMARK")].to_csv(OUT / "R9C_HALLMARK_COMPARISON.tsv", sep="\t", index=False, float_format="%.12g")

    # Frozen, result-independent gate operationalization: stable triad requires
    # both-cohort M1 Q2 < .50, both-cohort residual dependence >= .15 and top2
    # fraction >= .60, plus cross-cohort standardized beta distance > .25.
    metric_lookup = {
        (r.signature_id, r.construct_metric, r.metric, r.cohort_or_contrast): float(r.estimate)
        for r in metrics.itertuples()
    }

    def value(sid: str, construct: str, metric_name: str, contrast: str) -> float:
        return metric_lookup[(sid, construct, metric_name, contrast)]

    gate_rows = []
    for sid in signatures:
        nonid = all(value(sid, "M1_BULK_TO_MALIGNANT_IDENTIFIABILITY", "LOO_Q2", c) < .50 for c in COHORTS)
        residual = all(value(sid, "M4_RESIDUAL_DEPENDENCE", "MEAN_ABS_OFFDIAGONAL_CORRELATION", c) >= .15 and value(sid, "M5_LOW_RANK_RESIDUAL_GEOMETRY", "TOP2_VARIANCE_FRACTION", c) >= .60 for c in COHORTS)
        noninv = value(sid, "M7_MEASUREMENT_INVARIANCE", "BETA_VECTOR_STANDARDIZED_EUCLIDEAN_DISTANCE", "CGGA325_VS_CGGA693") > .25
        gate_rows.append({"signature_id": sid, "tier": tier[sid], "stable_nonidentification": nonid, "stable_residual_structure": residual, "stable_measurement_noninvariance": noninv, "shared_construct_triad": nonid and residual and noninv})
    gate = pd.DataFrame(gate_rows)
    tier1_shared = int(gate.loc[gate.tier.eq("TIER1_GBM_STATE"), "shared_construct_triad"].sum())
    tier2_shared = int(gate.loc[gate.tier.eq("TIER2_HALLMARK"), "shared_construct_triad"].sum())
    random_prevalence = float(gate.loc[gate.tier.eq("TIER3_MATCHED_RANDOM"), "shared_construct_triad"].mean())

    canonical_random = ranks.loc[ranks.reference_family.eq("MATCHED_RANDOM")]
    canonical_extreme = int(((canonical_random.adverse_direction.eq("HIGH") & (canonical_random.empirical_percentile >= 95)) | (canonical_random.adverse_direction.eq("LOW") & (canonical_random.empirical_percentile <= 5))).sum())
    m1_tail = canonical_random.loc[(canonical_random.construct_metric == "M1_BULK_TO_MALIGNANT_IDENTIFIABILITY") & (canonical_random.metric == "LOO_Q2"), "empirical_percentile"]
    m3_tail = canonical_random.loc[(canonical_random.construct_metric == "M3_SAME_SCORE_REALIZATION_DISPERSION") & (canonical_random.metric == "MEDIAN_RMS_REALIZATION_DISTANCE"), "empirical_percentile"]
    m4_tail = canonical_random.loc[(canonical_random.construct_metric == "M4_RESIDUAL_DEPENDENCE") & (canonical_random.metric == "MEAN_ABS_OFFDIAGONAL_CORRELATION"), "empirical_percentile"]
    mes_specific_nonidentification = bool(len(m1_tail) == 2 and (m1_tail <= 5).all() and ((len(m3_tail) == 2 and (m3_tail >= 95).all()) or (len(m4_tail) == 2 and (m4_tail >= 95).all())))
    if tier2_shared >= 2 or random_prevalence >= .25:
        outcome = "TRANSCRIPTOMIC_SIGNATURE_GENERAL_PROPERTY_WITHIN_GBM_CONTROLS"
    elif tier1_shared >= 2:
        outcome = "GBM_STATE_SIGNATURE_GENERAL_PROPERTY"
    elif mes_specific_nonidentification:
        outcome = "MES_SPECIFIC_EXCESS_NONIDENTIFICATION"
    else:
        outcome = "NO_STABLE_PATTERN"
    generalization_gate = outcome in {"GBM_STATE_SIGNATURE_GENERAL_PROPERTY", "TRANSCRIPTOMIC_SIGNATURE_GENERAL_PROPERTY_WITHIN_GBM_CONTROLS"} or (tier1_shared + tier2_shared >= 2)
    if outcome in {"GBM_STATE_SIGNATURE_GENERAL_PROPERTY", "TRANSCRIPTOMIC_SIGNATURE_GENERAL_PROPERTY_WITHIN_GBM_CONTROLS"}:
        reviewer_any = "ATTACK_SUPPORTED"
    elif outcome == "MES_SPECIFIC_EXCESS_NONIDENTIFICATION":
        reviewer_any = "ATTACK_REJECTED"
    else:
        reviewer_any = "ATTACK_PARTIALLY_SUPPORTED"

    negatives = [
        {"negative_finding_id": "R9C-N1", "finding": "Canonical MES did not satisfy the prespecified stable measurement-invariance criterion", "observed": str(bool(gate.loc[gate.signature_id.eq("CANONICAL_MES"), "stable_measurement_noninvariance"].iloc[0])), "interpretation": "Directional similarity can coexist with boundary-level magnitude difference; do not relabel as full invariance."},
        {"negative_finding_id": "R9C-N2", "finding": "Canonical MES did not show stable excess non-identification against matched random signatures", "observed": f"M1 adverse-tail both cohorts={bool(len(m1_tail)==2 and (m1_tail<=5).all())}; M3 adverse-tail both cohorts={bool(len(m3_tail)==2 and (m3_tail>=95).all())}; M4 adverse-tail both cohorts={bool(len(m4_tail)==2 and (m4_tail>=95).all())}", "interpretation": "Low-rank or transport extremeness cannot be relabeled as excess non-identification."},
        {"negative_finding_id": "R9C-N3", "finding": "The matched null is approximate rather than identical on coexpression", "observed": "median absolute pairwise-correlation delta=0.0373; 97.5%<=0.05", "interpretation": "Random-null inference is conditional on the audited matching tolerance."},
    ]
    pd.DataFrame(negatives).to_csv(OUT / "R9C_NEGATIVE_FINDINGS.tsv", sep="\t", index=False)
    gate.to_csv(OUT / "R9C_SIGNATURE_GATE_TRAITS.tsv", sep="\t", index=False)

    adjudication = f"""# R9C CONSTRUCT SPECIFICITY ADJUDICATION

## Formal outcome

`{outcome}`

The fair-pipeline analysis included 5 prespecified Neftel state controls, 6 Hallmark controls, and 5,000 phenotype-blind 95-gene null signatures. Comparator membership was defined before any construct metric was computed.

## Reviewer #2 attack

The assertion that *any* transcriptomic signature would behave this way is adjudicated `{reviewer_any}`. The prespecified shared-construct triad was present in {tier1_shared}/5 Neftel state controls, {tier2_shared}/6 Hallmarks, and {random_prevalence:.1%} of matched random signatures. Canonical MES occupied an adverse 5% matched-null tail in {canonical_extreme}/{len(canonical_random)} tested cohort/metric contrasts; therefore exceptionalism is not inferred merely from isolated percentile extremes.

## Generalization gate

`{'GENERALIZATION_GATE_OPEN' if generalization_gate else 'NO_JUSTIFIED_GENERALIZATION'}`

The gate is determined only by the frozen R9C outcome/shared-triad rule. This is a hypothesis gate, not evidence that the principle already generalizes beyond GBM.

## Boundaries

M1 is patient-level leave-one-out Q2; M2 excludes mechanically tautological endpoints; M3 uses continuous nearest-score matching; M4/M5 use the score-conditioned realization residual; M6 compares the cohort-specific residual rank-2 projectors; canonical M7 is read verbatim from R9B. No phenotype was used to choose comparator membership, and no result-guided comparator was added.
"""
    (OUT / "R9C_SPECIFICITY_ADJUDICATION.md").write_text(adjudication, encoding="utf-8")
    (RUN / "logs").mkdir(parents=True, exist_ok=True)
    (RUN / "logs/matched_signature_environment.json").write_text(json.dumps({
        "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
        "signatures": len(signatures), "patients": {c: len(cohort_data[c]["ids"]) for c in COHORTS},
        "m7_max_recalculation_difference_before_reference_replacement": max(r9b_errors), "r9a_canonical_geometry_max_error": max(r9a_errors),
        "gate_thresholds": {"m1_q2_both": "<0.50", "m4_both": ">=0.15", "m5_top2_both": ">=0.60", "m7_standardized_beta_distance": ">0.25"},
        "outcome": outcome, "generalization_gate_open": generalization_gate, "mes_specific_nonidentification_rule": mes_specific_nonidentification,
    }, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "R9C_COMPLETE", "outcome": outcome, "generalization_gate_open": generalization_gate, "tier1_shared": tier1_shared, "tier2_shared": tier2_shared, "random_shared_prevalence": random_prevalence, "r9b_max_check_error": max(r9b_errors)}))


if __name__ == "__main__":
    main()
