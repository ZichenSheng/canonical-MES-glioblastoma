#!/usr/bin/env python3
"""Run the fixed-rank, rotation-invariant V1.7.1 formal analysis."""

from __future__ import annotations

import itertools
import json
import platform
import sys
import os
from pathlib import Path

DATA_ROOT = Path(os.environ.get("GBM_MES_DATA_ROOT", "data"))
RESULT_ROOT = Path(os.environ.get("GBM_MES_OUTPUT_ROOT", "results"))
REPO_ROOT = Path(os.environ.get("GBM_MES_REPO_ROOT", Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


RUN = RESULT_ROOT / "rank2_subspace"
V17 = RESULT_ROOT / "structured_realization"
V16 = RESULT_ROOT / "patient_realization"
AXES = ["MALIGNANT_STATE", "MYELOID", "HYPOXIA", "MATRIX", "VASCULAR_STROMAL"]
COLS = [f"residual_z_{a}" for a in AXES]
N_BOOT = 5000
N_IND_PERM = 5000
N_GRASS = 100000
SEEDS = {
    "CGGA325_INTERNAL": 2026081711,
    "CGGA693_INTERNAL": 2026081712,
    "VALIDATION_ONLY": 2026081713,
    "TWO_COHORT": 2026081714,
    "INDEPENDENT_AXIS": 2026081715,
    "GRASSMANN": 2026081716,
}


def write(df: pd.DataFrame, rel: str) -> None:
    path = RUN / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="NA")


def pca_geometry(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xc = np.asarray(x, dtype=float) - np.mean(x, axis=0, keepdims=True)
    cov = (xc.T @ xc) / (xc.shape[0] - 1)
    values, vectors = np.linalg.eigh(cov)
    order = np.argsort(values)[::-1]
    values = values[order]
    vectors = vectors[:, order]
    return values, vectors, cov


def subspace_metrics(u: np.ndarray, v: np.ndarray) -> dict:
    singular = np.linalg.svd(u.T @ v, compute_uv=False)
    singular = np.clip(np.sort(singular)[::-1], 0.0, 1.0)
    angles = np.degrees(np.arccos(singular))
    sqsum = float(np.sum(singular**2))
    return {
        "theta1_degrees": float(angles[0]),
        "theta2_degrees": float(angles[1]),
        "cos_theta1": float(singular[0]),
        "cos_theta2": float(singular[1]),
        "S_proj": sqsum / 2.0,
        "D_chordal": float(np.sqrt(max(0.0, 2.0 - sqsum))),
    }


def bootstrap_subspaces(x: np.ndarray, reference: np.ndarray, seed: int, include_leverage: bool) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    n = x.shape[0]
    for rep in range(1, N_BOOT + 1):
        idx = rng.integers(0, n, n)
        _, vectors, _ = pca_geometry(x[idx])
        u = vectors[:, :2]
        row = {"replicate": rep, **subspace_metrics(reference, u)}
        if include_leverage:
            lev = np.diag(u @ u.T)
            row.update({f"leverage_{axis}": float(lev[j]) for j, axis in enumerate(AXES)})
        rows.append(row)
    return pd.DataFrame(rows)


def pct(x: pd.Series | np.ndarray, q: float) -> float:
    return float(np.quantile(np.asarray(x, dtype=float), q))


def corr_effects(x: np.ndarray, y: np.ndarray) -> dict:
    pr = pearsonr(x, y)
    sr = spearmanr(x, y)
    denom = np.linalg.norm(x) * np.linalg.norm(y)
    return {
        "pearson": float(pr.statistic),
        "pearson_p_descriptive": float(pr.pvalue),
        "spearman": float(sr.statistic),
        "spearman_p_descriptive": float(sr.pvalue),
        "cosine": float(np.dot(x, y) / denom),
    }


def main() -> None:
    validation = json.loads((RUN / "00_governance/freeze/validation_report.json").read_text())
    if validation["status"] != "PASS" or validation["blocker_count"] != 0:
        raise RuntimeError("Freeze validation did not pass")

    residual = pd.read_csv(V16 / "02_fingerprint/MES_CONDITIONED_FINGERPRINT.tsv", sep="\t")
    matrices = {}
    for cohort, n in [("CGGA325", 73), ("CGGA693", 108)]:
        frame = residual.loc[residual.cohort == cohort].copy()
        if len(frame) != n or frame.patient_id.nunique() != n or frame[COLS].isna().any().any():
            raise RuntimeError(f"{cohort} reference mismatch")
        matrices[cohort] = frame[COLS].to_numpy(float)

    # A-B: complete covariance geometry and identification diagnostics.
    eigens = []
    diagnostics = []
    refs = {}
    covs = {}
    for cohort, x in matrices.items():
        values, vectors, cov = pca_geometry(x)
        refs[cohort] = vectors[:, :2]
        covs[cohort] = cov
        fractions = values / values.sum()
        for j in range(5):
            eigens.append({
                "cohort": cohort, "component": f"PC{j+1}", "eigenvalue": values[j],
                "variance_fraction": fractions[j], "cumulative_variance_fraction": fractions[:j+1].sum(),
                **{f"loading_{axis}": vectors[k, j] for k, axis in enumerate(AXES)},
            })
        diagnostics.append({
            "cohort": cohort,
            "lambda1": values[0], "lambda2": values[1], "lambda3": values[2], "lambda4": values[3], "lambda5": values[4],
            "VF_PC1": fractions[0], "VF_PC2": fractions[1], "VF_PC3": fractions[2],
            "VF_PC1_PC2": fractions[:2].sum(),
            "G12": (values[0] - values[1]) / values[0],
            "G23": (values[1] - values[2]) / values[1],
            "trace_covariance": np.trace(cov),
            "rank2_orthonormality_max_abs_error": np.max(np.abs(refs[cohort].T @ refs[cohort] - np.eye(2))),
        })
    write(pd.DataFrame(eigens), "02_eigenspectrum/COHORT_EIGENSPECTRUM.tsv")
    write(pd.DataFrame(diagnostics), "02_eigenspectrum/RANK2_IDENTIFICATION_DIAGNOSTICS.tsv")

    # F-G: within-cohort patient bootstrap stability and leverage uncertainty.
    internal_tables = {}
    internal_summary = []
    leverage_boot_rows = []
    for cohort in ["CGGA325", "CGGA693"]:
        boot = bootstrap_subspaces(matrices[cohort], refs[cohort], SEEDS[f"{cohort}_INTERNAL"], True)
        internal_tables[cohort] = boot
        write(boot, f"03_within_cohort_stability/{cohort}_RANK2_BOOTSTRAP.tsv")
        low = pct(boot.S_proj, 0.025)
        median_s = pct(boot.S_proj, 0.5)
        median_t2 = pct(boot.theta2_degrees, 0.5)
        g1 = median_s >= 0.80
        g2 = low >= 0.60
        g3 = median_t2 <= 25.0
        internal_summary.append({
            "cohort": cohort,
            "median_S_proj": median_s,
            "S_proj_ci_low": low,
            "S_proj_ci_high": pct(boot.S_proj, 0.975),
            "median_theta1_degrees": pct(boot.theta1_degrees, 0.5),
            "median_theta2_degrees": median_t2,
            "theta2_P95_degrees": pct(boot.theta2_degrees, 0.95),
            "proportion_S_proj_ge_0_80": float(np.mean(boot.S_proj >= 0.80)),
            "proportion_theta2_le_25deg": float(np.mean(boot.theta2_degrees <= 25.0)),
            "gate_G1_median_S_ge_0_80": g1,
            "gate_G2_lower_CI_ge_0_60": g2,
            "gate_G3_median_theta2_le_25": g3,
            "internal_stability_status": "RANK2_INTERNAL_STABILITY_PASS" if g1 and g2 and g3 else "RANK2_INTERNAL_STABILITY_FAIL",
            "bootstrap_replicates": N_BOOT,
        })
        for _, row in boot.iterrows():
            for axis in AXES:
                leverage_boot_rows.append({"cohort": cohort, "replicate": int(row.replicate), "axis": axis, "leverage": row[f"leverage_{axis}"]})
    internal_summary_df = pd.DataFrame(internal_summary)
    write(internal_summary_df, "03_within_cohort_stability/RANK2_INTERNAL_STABILITY_SUMMARY.tsv")
    leverage_boot = pd.DataFrame(leverage_boot_rows)
    write(leverage_boot, "06_rotation_invariant_interpretation/SUBSPACE_LEVERAGE_BOOTSTRAP.tsv")

    # H: observed fixed discovery-to-validation comparison.
    u325, u693 = refs["CGGA325"], refs["CGGA693"]
    observed = subspace_metrics(u325, u693)
    angle_rows = [
        {"comparison": "CGGA325_DISCOVERY_VS_CGGA693_VALIDATION", "principal_angle": "theta1", "degrees": observed["theta1_degrees"], "cosine": observed["cos_theta1"]},
        {"comparison": "CGGA325_DISCOVERY_VS_CGGA693_VALIDATION", "principal_angle": "theta2", "degrees": observed["theta2_degrees"], "cosine": observed["cos_theta2"]},
    ]
    write(pd.DataFrame(angle_rows), "04_cross_cohort_subspace/PRINCIPAL_ANGLES.tsv")
    p325, p693 = u325 @ u325.T, u693 @ u693.T
    projection_qa = {
        "comparison": "CGGA325_DISCOVERY_VS_CGGA693_VALIDATION", **observed,
        "trace_P325_P693_over_2_direct": float(np.trace(p325 @ p693) / 2.0),
        "P325_symmetry_max_error": float(np.max(np.abs(p325 - p325.T))),
        "P693_symmetry_max_error": float(np.max(np.abs(p693 - p693.T))),
        "P325_idempotence_max_error": float(np.max(np.abs(p325 @ p325 - p325))),
        "P693_idempotence_max_error": float(np.max(np.abs(p693 @ p693 - p693))),
        "rank_P325": int(np.linalg.matrix_rank(p325, tol=1e-10)),
        "rank_P693": int(np.linalg.matrix_rank(p693, tol=1e-10)),
    }
    write(pd.DataFrame([projection_qa]), "04_cross_cohort_subspace/PROJECTION_SIMILARITY.tsv")

    # I: validation-only primary uncertainty (U325 fixed).
    validation_boot = bootstrap_subspaces(matrices["CGGA693"], u325, SEEDS["VALIDATION_ONLY"], False)
    write(validation_boot, "04_cross_cohort_subspace/VALIDATION_BOOTSTRAP.tsv")

    # J: secondary two-cohort patient bootstrap sensitivity.
    rng = np.random.default_rng(SEEDS["TWO_COHORT"])
    x325, x693 = matrices["CGGA325"], matrices["CGGA693"]
    two_rows = []
    for rep_i in range(1, N_BOOT + 1):
        _, b325, _ = pca_geometry(x325[rng.integers(0, len(x325), len(x325))])
        _, b693, _ = pca_geometry(x693[rng.integers(0, len(x693), len(x693))])
        two_rows.append({"replicate": rep_i, **subspace_metrics(b325[:, :2], b693[:, :2])})
    two_boot = pd.DataFrame(two_rows)
    write(two_boot, "04_cross_cohort_subspace/TWO_COHORT_BOOTSTRAP.tsv")

    # K: exact 5! coordinate-label permutation null.
    exact_rows = []
    identity = tuple(range(5))
    for perm in itertools.permutations(range(5)):
        metrics = subspace_metrics(u325, u693[list(perm), :])
        exact_rows.append({
            "permutation_index": len(exact_rows) + 1,
            "permutation_zero_based": ";".join(map(str, perm)),
            "permuted_axis_labels_in_discovery_coordinate_order": ";".join(AXES[j] for j in perm),
            "identity_permutation": perm == identity,
            **metrics,
        })
    exact = pd.DataFrame(exact_rows)
    write(exact, "05_null_models/EXACT_AXIS_LABEL_PERMUTATION.tsv")

    # L: independent-axis patient-label permutation null in CGGA693.
    rng = np.random.default_rng(SEEDS["INDEPENDENT_AXIS"])
    independent_rows = []
    for rep_i in range(1, N_IND_PERM + 1):
        xp = np.column_stack([x693[rng.permutation(len(x693)), j] for j in range(5)])
        _, vp, _ = pca_geometry(xp)
        independent_rows.append({"replicate": rep_i, **subspace_metrics(u325, vp[:, :2])})
    independent = pd.DataFrame(independent_rows)
    write(independent, "05_null_models/INDEPENDENT_AXIS_PERMUTATION.tsv")

    # M: random Grassmann geometry sensitivity in batches.
    rng = np.random.default_rng(SEEDS["GRASSMANN"])
    grass_rows = []
    generated = 0
    batch_size = 5000
    while generated < N_GRASS:
        batch = min(batch_size, N_GRASS - generated)
        gaussian = rng.normal(size=(batch, 5, 2))
        q, _ = np.linalg.qr(gaussian)
        for j in range(batch):
            grass_rows.append({"replicate": generated + j + 1, **subspace_metrics(u325, q[j, :, :2])})
        generated += batch
    grass = pd.DataFrame(grass_rows)
    write(grass, "05_null_models/RANDOM_GRASSMANN_NULL.tsv")

    # Null and bootstrap summaries.
    exact_p = float(np.mean(exact.S_proj >= observed["S_proj"] - 1e-15))
    independent_p = float((1 + np.sum(independent.S_proj >= observed["S_proj"] - 1e-15)) / (N_IND_PERM + 1))
    grass_p = float(np.mean(grass.S_proj >= observed["S_proj"] - 1e-15))
    null_summary = pd.DataFrame([
        {"null_model": "EXACT_AXIS_LABEL_PERMUTATION", "iterations": 120, "observed_S_proj": observed["S_proj"], "null_mean": exact.S_proj.mean(), "null_median": exact.S_proj.median(), "null_P95": pct(exact.S_proj, .95), "null_P99": pct(exact.S_proj, .99), "max_nonidentity_S_proj": exact.loc[~exact.identity_permutation, "S_proj"].max(), "tail_probability": exact_p, "p_definition": "complete enumeration; count>=observed / 120"},
        {"null_model": "INDEPENDENT_AXIS_PATIENT_PERMUTATION", "iterations": N_IND_PERM, "observed_S_proj": observed["S_proj"], "null_mean": independent.S_proj.mean(), "null_median": independent.S_proj.median(), "null_P95": pct(independent.S_proj, .95), "null_P99": pct(independent.S_proj, .99), "max_nonidentity_S_proj": np.nan, "tail_probability": independent_p, "p_definition": "finite (1+count>=observed)/5001"},
        {"null_model": "RANDOM_GRASSMANN_GEOMETRY", "iterations": N_GRASS, "observed_S_proj": observed["S_proj"], "null_mean": grass.S_proj.mean(), "null_median": grass.S_proj.median(), "null_P95": pct(grass.S_proj, .95), "null_P99": pct(grass.S_proj, .99), "max_nonidentity_S_proj": np.nan, "tail_probability": grass_p, "p_definition": "empirical geometric tail count/100000"},
    ])
    write(null_summary, "05_null_models/NULL_MODEL_SUMMARY.tsv")

    # N-O: formal leverage and off-diagonal projector geometry.
    leverage_rows = []
    for cohort, p in [("CGGA325", p325), ("CGGA693", p693)]:
        cohort_boot = leverage_boot.loc[leverage_boot.cohort == cohort]
        for j, axis in enumerate(AXES):
            b = cohort_boot.loc[cohort_boot.axis == axis, "leverage"]
            leverage_rows.append({
                "cohort": cohort, "axis": axis, "formal_leverage": p[j, j],
                "bootstrap_median": b.median(), "bootstrap_ci_low": pct(b, .025), "bootstrap_ci_high": pct(b, .975),
                "formal_leverage_sum_within_cohort": float(np.trace(p)),
            })
    leverage = pd.DataFrame(leverage_rows)
    write(leverage, "06_rotation_invariant_interpretation/SUBSPACE_LEVERAGE.tsv")

    v17_corr = pd.read_csv(V17 / "02_residual_structure/RESIDUAL_SPEARMAN_CORRELATIONS.tsv", sep="\t")
    off_rows = []
    for cohort, p in [("CGGA325", p325), ("CGGA693", p693)]:
        for a in range(5):
            for b in range(a + 1, 5):
                frozen = v17_corr.loc[(v17_corr.cohort == cohort) & (v17_corr.axis_a == AXES[a]) & (v17_corr.axis_b == AXES[b])]
                if len(frozen) != 1:
                    raise RuntimeError(f"V1.7 correlation pair mismatch: {cohort} {AXES[a]} {AXES[b]}")
                off_rows.append({
                    "cohort": cohort, "axis_a": AXES[a], "axis_b": AXES[b],
                    "projection_offdiagonal": p[a, b], "v1_7_frozen_spearman_rho": frozen.iloc[0].spearman_rho,
                    "interpretation_boundary": "ROTATION_INVARIANT_GEOMETRY_NOT_CAUSAL_INTERACTION",
                })
    off = pd.DataFrame(off_rows)
    write(off, "06_rotation_invariant_interpretation/PROJECTION_OFFDIAGONAL_GEOMETRY.tsv")

    l325 = leverage.loc[leverage.cohort == "CGGA325", "formal_leverage"].to_numpy()
    l693 = leverage.loc[leverage.cohort == "CGGA693", "formal_leverage"].to_numpy()
    o325 = off.loc[off.cohort == "CGGA325", "projection_offdiagonal"].to_numpy()
    o693 = off.loc[off.cohort == "CGGA693", "projection_offdiagonal"].to_numpy()
    concord_rows = []
    le = corr_effects(l325, l693)
    concord_rows.append({"object": "FIVE_AXIS_LEVERAGE_VECTOR", "n_elements": 5, **le, "mean_absolute_difference": float(np.mean(np.abs(l325 - l693))), "interpretation_boundary": "EFFECT_SIZE_ONLY_N5"})
    oe = corr_effects(o325, o693)
    concord_rows.append({"object": "TEN_OFFDIAGONAL_PROJECTION_ELEMENTS", "n_elements": 10, **oe, "mean_absolute_difference": float(np.mean(np.abs(o325 - o693))), "interpretation_boundary": "ROTATION_INVARIANT_GEOMETRY_NOT_CAUSAL_INTERACTION"})
    for cohort in ["CGGA325", "CGGA693"]:
        sub = off.loc[off.cohort == cohort]
        ce = corr_effects(sub.projection_offdiagonal.to_numpy(), sub.v1_7_frozen_spearman_rho.to_numpy())
        concord_rows.append({"object": f"PROJECTION_GEOMETRY_VS_V1_7_RHO_{cohort}", "n_elements": 10, **ce, "mean_absolute_difference": float(np.mean(np.abs(sub.projection_offdiagonal - sub.v1_7_frozen_spearman_rho))), "interpretation_boundary": "DESCRIPTIVE_ONLY; V1_7_CLASSIFICATION_UNCHANGED"})
    concord_df = pd.DataFrame(concord_rows)
    write(concord_df, "06_rotation_invariant_interpretation/CROSS_COHORT_ROTATION_INVARIANT_CONCORDANCE.tsv")

    # Q-R: explicit single-axis versus fixed rank-2 comparison and full gate.
    v17_pca325 = pd.read_csv(V17 / "03_realization_axis/CGGA325_DISCOVERY_PCA.tsv", sep="\t")
    v17_pca693 = pd.read_csv(V17 / "04_cross_cohort_validation/CGGA693_INDEPENDENT_PCA.tsv", sep="\t")
    v17_single = pd.read_csv(V17 / "04_cross_cohort_validation/LOADING_CONCORDANCE.tsv", sep="\t").iloc[0]
    v17_gate = pd.read_csv(V17 / "04_cross_cohort_validation/CROSS_COHORT_REALIZATION_GATE.tsv", sep="\t").iloc[0]
    val_ci_low = pct(validation_boot.S_proj, .025)
    all_internal = internal_summary_df.internal_stability_status.eq("RANK2_INTERNAL_STABILITY_PASS").all()
    gate_conditions = {
        "R1_both_internal_stability_pass": bool(all_internal),
        "R2_observed_S_proj_ge_0_80": observed["S_proj"] >= 0.80,
        "R3_theta2_le_25deg": observed["theta2_degrees"] <= 25.0,
        "R4_exact_axis_p_le_0_05": exact_p <= 0.05,
        "R5_independent_axis_p_le_0_01": independent_p <= 0.01,
        "R6_validation_bootstrap_lower_S_ge_0_65": val_ci_low >= 0.65,
    }
    gate_pass = all(gate_conditions.values())
    if not all_internal:
        verdict = "RANK2_REALIZATION_SUBSPACE_NOT_STABLY_IDENTIFIED"
    elif gate_pass:
        verdict = "STRUCTURED_MULTIDIMENSIONAL_REALIZATION_WITH_REPRODUCIBLE_RANK2_SUBSPACE"
    else:
        verdict = "STRUCTURED_MULTIDIMENSIONAL_REALIZATION_WITHOUT_CROSS_COHORT_SUBSPACE_REPLICATION"
    comparison = pd.DataFrame([
        {"model": "V1_7_SINGLE_PC1", "metric": "PC1_cosine_similarity", "value": v17_single.cosine_similarity, "status": "V1_7_SINGLE_AXIS_TRANSFERABILITY_FAIL"},
        {"model": "V1_7_SINGLE_PC1", "metric": "CGGA325_PC1_variance_fraction", "value": v17_pca325.loc[v17_pca325.component == "PC1", "variance_fraction"].iloc[0], "status": v17_gate.cross_cohort_gate_status},
        {"model": "V1_7_SINGLE_PC1", "metric": "CGGA325_PC2_variance_fraction", "value": v17_pca325.loc[v17_pca325.component == "PC2", "variance_fraction"].iloc[0], "status": v17_gate.cross_cohort_gate_status},
        {"model": "V1_7_1_RANK2_SUBSPACE", "metric": "CGGA325_rank2_cumulative_variance", "value": pd.DataFrame(diagnostics).loc[0, "VF_PC1_PC2"], "status": verdict},
        {"model": "V1_7_1_RANK2_SUBSPACE", "metric": "CGGA693_rank2_cumulative_variance", "value": pd.DataFrame(diagnostics).loc[1, "VF_PC1_PC2"], "status": verdict},
        {"model": "V1_7_1_RANK2_SUBSPACE", "metric": "theta1_degrees", "value": observed["theta1_degrees"], "status": verdict},
        {"model": "V1_7_1_RANK2_SUBSPACE", "metric": "theta2_degrees", "value": observed["theta2_degrees"], "status": verdict},
        {"model": "V1_7_1_RANK2_SUBSPACE", "metric": "S_proj", "value": observed["S_proj"], "status": verdict},
        {"model": "V1_7_1_RANK2_SUBSPACE", "metric": "exact_axis_label_p", "value": exact_p, "status": verdict},
        {"model": "V1_7_1_RANK2_SUBSPACE", "metric": "independent_axis_p", "value": independent_p, "status": verdict},
    ])
    write(comparison, "07_v17_integration/SINGLE_AXIS_VS_SUBSPACE_SUMMARY.tsv")

    gate_row = {
        **gate_conditions,
        "validation_only_S_proj_ci_low": val_ci_low,
        "subspace_validation_status": "SUBSPACE_VALIDATION_PASS" if gate_pass else "SUBSPACE_VALIDATION_FAIL",
        "scientific_verdict": verdict,
        "branch_stop": "PATIENT_REALIZATION_STATISTICAL_BRANCH_FROZEN_NO_V1_7_2",
    }
    write(pd.DataFrame([gate_row]), "07_v17_integration/SUBSPACE_VALIDATION_GATE.tsv")

    # Compact formal result used by finalizer.
    final = {
        "comparison": "CGGA325_DISCOVERY_VS_CGGA693_VALIDATION",
        **observed,
        "CGGA325_rank2_variance_fraction": pd.DataFrame(diagnostics).loc[0, "VF_PC1_PC2"],
        "CGGA693_rank2_variance_fraction": pd.DataFrame(diagnostics).loc[1, "VF_PC1_PC2"],
        "validation_bootstrap_median_S_proj": validation_boot.S_proj.median(),
        "validation_bootstrap_S_proj_ci_low": pct(validation_boot.S_proj, .025),
        "validation_bootstrap_S_proj_ci_high": pct(validation_boot.S_proj, .975),
        "validation_bootstrap_theta2_ci_low": pct(validation_boot.theta2_degrees, .025),
        "validation_bootstrap_theta2_ci_high": pct(validation_boot.theta2_degrees, .975),
        "two_cohort_bootstrap_median_S_proj": two_boot.S_proj.median(),
        "two_cohort_bootstrap_S_proj_ci_low": pct(two_boot.S_proj, .025),
        "two_cohort_bootstrap_S_proj_ci_high": pct(two_boot.S_proj, .975),
        "exact_axis_label_p": exact_p,
        "independent_axis_permutation_p": independent_p,
        "grassmann_geometric_tail_probability": grass_p,
        "subspace_validation_status": gate_row["subspace_validation_status"],
        "scientific_verdict": verdict,
    }
    write(pd.DataFrame([final]), "FINAL_REALIZATION_SUBSPACE_V1_7_1.tsv")

    session = {
        "python": sys.version, "platform": platform.platform(), "numpy": np.__version__, "pandas": pd.__version__,
        "scipy_available": True, "rank": 2, "axes": AXES, "bootstrap_iterations": N_BOOT,
        "independent_axis_permutations": N_IND_PERM, "grassmann_simulations": N_GRASS, "seeds": SEEDS,
    }
    (RUN / "session_info/FORMAL_ANALYSIS_SESSION_INFO.json").write_text(json.dumps(session, indent=2) + "\n")
    (RUN / "logs/formal_subspace_analysis.log").write_text(json.dumps({"status": "COMPLETE", "observed": observed, "gate": gate_row}, indent=2) + "\n")
    print(json.dumps({"observed": observed, "internal": internal_summary, "nulls": null_summary.to_dict(orient="records"), "gate": gate_row}, indent=2))


if __name__ == "__main__":
    main()
