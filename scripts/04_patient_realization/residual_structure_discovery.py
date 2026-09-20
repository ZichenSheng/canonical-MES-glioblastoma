#!/usr/bin/env python3
"""Formal Modules A-G: residual structure and locked CGGA325 PC1 discovery."""

from __future__ import annotations

import hashlib
import json
import math
import platform
from datetime import datetime, timezone
from itertools import combinations
import os
from pathlib import Path

DATA_ROOT = Path(os.environ.get("GBM_MES_DATA_ROOT", "data"))
RESULT_ROOT = Path(os.environ.get("GBM_MES_OUTPUT_ROOT", "results"))
REPO_ROOT = Path(os.environ.get("GBM_MES_REPO_ROOT", Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd
import scipy
from scipy.stats import pearsonr, spearmanr


RUN = RESULT_ROOT / "structured_realization"
V16 = RESULT_ROOT / "patient_realization"
AXES = ["MALIGNANT_STATE", "MYELOID", "HYPOXIA", "MATRIX", "VASCULAR_STROMAL"]
COLS = [f"residual_z_{axis}" for axis in AXES]
EXPECTED_N = {"CGGA325": 73, "CGGA693": 108}
B = 5000


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(frame: pd.DataFrame, relative: str) -> None:
    path = RUN / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA", lineterminator="\n")


def bh_adjust(values: np.ndarray) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out = np.empty_like(adjusted)
    out[order] = adjusted
    return out


def corr_matrix_spearman(x: np.ndarray) -> np.ndarray:
    result = spearmanr(x, axis=0).statistic
    return np.asarray(result, dtype=float)


def verify_manifest(run_dir: Path, manifest: Path) -> tuple[int, list[str]]:
    mismatches = []
    checked = 0
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected, relative = line.split("  ", 1)
        target = run_dir / relative
        checked += 1
        if not target.is_file() or sha256(target) != expected:
            mismatches.append(relative)
    return checked, mismatches


def pca(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    centered = x - x.mean(axis=0, keepdims=True)
    _, singular, vt = np.linalg.svd(centered, full_matrices=False)
    eigenvalues = singular**2 / (x.shape[0] - 1)
    loadings = vt.T
    scores = centered @ loadings
    return eigenvalues, loadings, scores


def orient_primary(loadings: np.ndarray) -> tuple[np.ndarray, str, str]:
    vector = loadings.copy()
    if abs(vector[0]) >= 0.10:
        anchor_index = 0
        rule = "MALIGNANT_STATE_POSITIVE"
    else:
        anchor_index = int(np.argmax(np.abs(vector)))
        rule = "ABSOLUTE_LARGEST_LOADING_POSITIVE_FALLBACK"
    if vector[anchor_index] < 0:
        vector *= -1
    return vector, AXES[anchor_index], rule


def main() -> None:
    started = datetime.now(timezone.utc)
    seeds = json.loads((RUN / "00_governance/V1_7_RANDOM_SEEDS.json").read_text())["seeds"]
    protocol_sha = sha256(RUN / "00_governance/GBM_MES_STRUCTURED_REALIZATION_PROTOCOL_V1_7_FROZEN.md")
    validation = json.loads((RUN / "00_governance/freeze/validation_report.json").read_text())
    if validation.get("status") != "PASS" or validation.get("blocker_count") != 0:
        raise RuntimeError("Formal analysis blocked: freeze validator did not pass")

    checked, failed = verify_manifest(V16, V16 / "FINAL_CHECKSUM_MANIFEST_V1_6.sha256")
    if failed:
        raise RuntimeError(f"V1.6 checksum failures: {failed[:5]}")
    write_tsv(pd.DataFrame([{"run": str(V16), "files_checked": checked, "files_failed": len(failed), "status": "PASS_READ_ONLY_INTEGRITY"}]), "01_reference/V1_6_FULL_CHECKSUM_RECHECK.tsv")

    residual = pd.read_csv(V16 / "02_fingerprint/MES_CONDITIONED_FINGERPRINT.tsv", sep="\t")
    master = pd.read_csv(V16 / "01_reference/PATIENT_REALIZATION_MASTER.tsv", sep="\t")
    qa_rows = []
    cohort_data: dict[str, pd.DataFrame] = {}
    for cohort, expected in EXPECTED_N.items():
        data = residual.loc[residual.cohort == cohort].copy().sort_values("patient_id").reset_index(drop=True)
        auth = pd.read_csv(RUN / f"00_governance/freeze/membership_{cohort.lower()}_patients_v1_7.tsv", sep="\t", dtype=str)
        exact_membership = set(data.patient_id.astype(str)) == set(auth.patient_id.astype(str))
        checks = {
            "n": len(data),
            "expected_n": expected,
            "duplicate_patient_ids": int(data.patient_id.duplicated().sum()),
            "missing_residual_values": int(data[COLS].isna().sum().sum()),
            "max_abs_axis_mean": float(data[COLS].mean().abs().max()),
            "max_abs_axis_sd_minus_1": float((data[COLS].std(ddof=1) - 1).abs().max()),
            "exact_membership": exact_membership,
            "cross_fitted_provenance": "PASS_LINEAR_LOPO_COLUMNS_AND_SOURCE_CODE",
        }
        status = "PASS" if len(data) == expected and checks["duplicate_patient_ids"] == 0 and checks["missing_residual_values"] == 0 and checks["max_abs_axis_mean"] < 1e-10 and checks["max_abs_axis_sd_minus_1"] < 1e-10 and exact_membership else "FAIL"
        qa_rows.append({"cohort": cohort, **checks, "status": status})
        cohort_data[cohort] = data
    input_qa = pd.DataFrame(qa_rows)
    write_tsv(input_qa, "01_reference/RESIDUAL_FINGERPRINT_INPUT_QA.tsv")
    if not input_qa.status.eq("PASS").all():
        raise RuntimeError("STOP_REFERENCE_OR_FINGERPRINT_MISMATCH")

    # Modules A-C: all ten primary Spearman relationships and Pearson sensitivity.
    correlation_rows = []
    sensitivity_rows = []
    bootstrap_rows = []
    bootstrap_store: dict[tuple[str, str, str], np.ndarray] = {}
    for cohort_index, (cohort, data) in enumerate(cohort_data.items()):
        x = data[COLS].to_numpy(float)
        pairs = list(combinations(range(5), 2))
        point_p = []
        point_values = []
        pearson_p = []
        pearson_values = []
        for a, b in pairs:
            rho, p = spearmanr(x[:, a], x[:, b])
            r, pp = pearsonr(x[:, a], x[:, b])
            point_values.append(float(rho)); point_p.append(float(p))
            pearson_values.append(float(r)); pearson_p.append(float(pp))

        rng = np.random.default_rng(seeds["residual_correlation_bootstrap"] + cohort_index * 10000)
        boot_cube = np.empty((B, len(pairs)), dtype=float)
        for rep in range(B):
            idx = rng.integers(0, len(x), len(x))
            cm = corr_matrix_spearman(x[idx, :])
            for pi, (a, b) in enumerate(pairs):
                boot_cube[rep, pi] = cm[a, b]
        fdr = bh_adjust(np.array(point_p))
        pearson_fdr = bh_adjust(np.array(pearson_p))
        for pi, (a, b) in enumerate(pairs):
            axis_a, axis_b = AXES[a], AXES[b]
            boot = boot_cube[:, pi]
            bootstrap_store[(cohort, axis_a, axis_b)] = boot
            correlation_rows.append({
                "cohort": cohort, "axis_a": axis_a, "axis_b": axis_b, "n": len(x),
                "spearman_rho": point_values[pi], "bootstrap_ci_low": np.quantile(boot, 0.025), "bootstrap_ci_high": np.quantile(boot, 0.975),
                "nominal_p": point_p[pi], "BH_FDR": fdr[pi], "sign": "POSITIVE" if point_values[pi] > 0 else "NEGATIVE" if point_values[pi] < 0 else "ZERO",
                "bootstrap_replicates": B,
            })
            sensitivity_rows.append({"cohort": cohort, "axis_a": axis_a, "axis_b": axis_b, "method": "PEARSON", "estimate": pearson_values[pi], "nominal_p": pearson_p[pi], "BH_FDR": pearson_fdr[pi], "status": "EVALUATED"})
            sensitivity_rows.append({"cohort": cohort, "axis_a": axis_a, "axis_b": axis_b, "method": "BIWEIGHT_MIDCORRELATION", "estimate": np.nan, "nominal_p": np.nan, "BH_FDR": np.nan, "status": "NOT_RUN_DEPENDENCY_NOT_PRECERTIFIED"})
            for rep, value in enumerate(boot, start=1):
                bootstrap_rows.append({"cohort": cohort, "axis_a": axis_a, "axis_b": axis_b, "replicate": rep, "spearman_rho": value})

    correlations = pd.DataFrame(correlation_rows)
    write_tsv(correlations, "02_residual_structure/RESIDUAL_SPEARMAN_CORRELATIONS.tsv")
    write_tsv(pd.DataFrame(bootstrap_rows), "02_residual_structure/RESIDUAL_CORRELATION_BOOTSTRAP.tsv")
    write_tsv(pd.DataFrame(sensitivity_rows), "02_residual_structure/RESIDUAL_CORRELATION_SENSITIVITY.tsv")

    # Cross-cohort direction, strong replication, and descriptive Fisher-z pooling.
    replication_rows = []
    for a, b in combinations(AXES, 2):
        r325 = correlations.loc[(correlations.cohort == "CGGA325") & (correlations.axis_a == a) & (correlations.axis_b == b)].iloc[0]
        r693 = correlations.loc[(correlations.cohort == "CGGA693") & (correlations.axis_a == a) & (correlations.axis_b == b)].iloc[0]
        same_sign = np.sign(r325.spearman_rho) == np.sign(r693.spearman_rho) and r325.spearman_rho != 0 and r693.spearman_rho != 0
        effect_gate = abs(r325.spearman_rho) >= 0.20 and abs(r693.spearman_rho) >= 0.20
        direction_replicated = bool(same_sign and effect_gate)
        if same_sign and r325.spearman_rho > 0:
            other_not_opposite = r325.bootstrap_ci_high >= 0 and r693.bootstrap_ci_high >= 0
            classification = "CO_ENRICHMENT"
        elif same_sign and r325.spearman_rho < 0:
            other_not_opposite = r325.bootstrap_ci_low <= 0 and r693.bootstrap_ci_low <= 0
            classification = "COMPENSATORY_OR_SUBSTITUTIVE_PATTERN"
        else:
            other_not_opposite = False
            classification = "NOT_DIRECTION_REPLICATED"
        strong = direction_replicated and min(r325.BH_FDR, r693.BH_FDR) < 0.05 and other_not_opposite

        rho_values = np.clip(np.array([r325.spearman_rho, r693.spearman_rho], float), -0.999999, 0.999999)
        weights = np.array([r325.n - 3, r693.n - 3], float)
        zs = np.arctanh(rho_values)
        pooled_z = np.sum(weights * zs) / np.sum(weights)
        se = math.sqrt(1 / np.sum(weights))
        pooled_rho = math.tanh(pooled_z)
        pooled_low = math.tanh(pooled_z - 1.96 * se)
        pooled_high = math.tanh(pooled_z + 1.96 * se)
        q = float(np.sum(weights * (zs - pooled_z) ** 2))
        i2 = max(0.0, (q - 1) / q) * 100 if q > 0 else 0.0
        replication_rows.append({
            "axis_a": a, "axis_b": b, "rho_CGGA325": r325.spearman_rho, "ci325_low": r325.bootstrap_ci_low, "ci325_high": r325.bootstrap_ci_high, "FDR_CGGA325": r325.BH_FDR,
            "rho_CGGA693": r693.spearman_rho, "ci693_low": r693.bootstrap_ci_low, "ci693_high": r693.bootstrap_ci_high, "FDR_CGGA693": r693.BH_FDR,
            "direction_concordant": bool(same_sign), "effect_size_gate_both_abs_ge_0_20": bool(effect_gate), "direction_replication_status": "DIRECTION_REPLICATED" if direction_replicated else "NOT_DIRECTION_REPLICATED",
            "strong_replication_status": "STRONG_REPLICATED_STRUCTURE" if strong else "NOT_STRONG_REPLICATED_STRUCTURE", "pattern_class": classification if direction_replicated else "NOT_CLASSIFIED",
            "pooled_rho_fixed_descriptive": pooled_rho, "pooled_ci_low": pooled_low, "pooled_ci_high": pooled_high, "heterogeneity_Q_df1": q, "I2_percent_diagnostic": i2,
            "interpretation_boundary": "CORRELATION_DOES_NOT_ESTABLISH_MECHANISTIC_COMPENSATION",
        })
    replication = pd.DataFrame(replication_rows)
    write_tsv(replication, "02_residual_structure/RESIDUAL_STRUCTURE_REPLICATION.tsv")

    # Module D: independent-axis label permutations.
    global_perm_rows = []
    global_summary_rows = []
    for cohort_index, (cohort, data) in enumerate(cohort_data.items()):
        x = data[COLS].to_numpy(float)
        cm = corr_matrix_spearman(x)
        observed = sum(cm[a, b] ** 2 for a, b in combinations(range(5), 2))
        rng = np.random.default_rng(seeds["global_dependence_permutation"] + cohort_index * 10000)
        null = np.empty(B, float)
        permuted = np.empty_like(x)
        for rep in range(B):
            for column in range(5):
                permuted[:, column] = x[rng.permutation(len(x)), column]
            pcm = corr_matrix_spearman(permuted)
            null[rep] = sum(pcm[a, b] ** 2 for a, b in combinations(range(5), 2))
        p = (1 + int(np.sum(null >= observed))) / (B + 1)
        global_summary_rows.append({"cohort": cohort, "n": len(x), "T_observed": observed, "null_median": np.median(null), "null_P95": np.quantile(null, 0.95), "finite_permutation_p": p, "permutations": B, "status": "STRUCTURED_DEPENDENCE_SUPPORTED" if p < 0.05 else "GLOBAL_DEPENDENCE_NOT_SUPPORTED"})
        for rep, value in enumerate(null, start=1):
            global_perm_rows.append({"cohort": cohort, "replicate": rep, "T_permuted": value})
    global_summary = pd.DataFrame(global_summary_rows)
    write_tsv(pd.DataFrame(global_perm_rows), "02_residual_structure/GLOBAL_DEPENDENCE_PERMUTATION.tsv")
    write_tsv(global_summary, "02_residual_structure/GLOBAL_DEPENDENCE_SUMMARY.tsv")
    level1 = global_summary.finite_permutation_p.lt(0.05).all() and replication.strong_replication_status.eq("STRONG_REPLICATED_STRUCTURE").any()
    write_tsv(pd.DataFrame([{"global_dependence_both_cohorts_p_lt_0_05": global_summary.finite_permutation_p.lt(0.05).all(), "n_strong_replicated_pairs": int(replication.strong_replication_status.eq("STRONG_REPLICATED_STRUCTURE").sum()), "level1_structured_residual_gate": "PASS" if level1 else "FAIL"}]), "02_residual_structure/RESIDUAL_STRUCTURE_LEVEL1_GATE.tsv")

    # Modules E-G: CGGA325 discovery PCA, PC1 only primary.
    discovery = cohort_data["CGGA325"]
    x325 = discovery[COLS].to_numpy(float)
    eigenvalues, loadings_all, scores_all = pca(x325)
    formal_pc1, anchor_axis, orientation_rule = orient_primary(loadings_all[:, 0])
    if np.dot(formal_pc1, loadings_all[:, 0]) < 0:
        loadings_all[:, 0] *= -1
        scores_all[:, 0] *= -1
    variance_fraction = eigenvalues / eigenvalues.sum()
    scree = pd.DataFrame({"component": [f"PC{i}" for i in range(1, 6)], "eigenvalue": eigenvalues, "variance_fraction": variance_fraction, "cumulative_variance_fraction": np.cumsum(variance_fraction), "primary_status": ["PRIMARY_PC1"] + ["DESCRIPTIVE_QA_ONLY"] * 4})
    write_tsv(scree, "03_realization_axis/CGGA325_DISCOVERY_PCA.tsv")
    write_tsv(scree.copy(), "03_realization_axis/CGGA325_PCA_SCREE_TABLE.tsv")
    all_loading_rows = []
    for ai, axis in enumerate(AXES):
        for pc in range(5):
            all_loading_rows.append({"axis": axis, "component": f"PC{pc+1}", "loading": loadings_all[ai, pc], "primary_status": "PRIMARY_PC1" if pc == 0 else "DESCRIPTIVE_QA_ONLY"})
    write_tsv(pd.DataFrame(all_loading_rows), "03_realization_axis/CGGA325_ALL_PCA_LOADINGS_QA.tsv")

    rng = np.random.default_rng(seeds["pca_loading_bootstrap"])
    boot_loadings = np.empty((B, 5), float)
    for rep in range(B):
        idx = rng.integers(0, len(x325), len(x325))
        _, lb, _ = pca(x325[idx, :])
        vector = lb[:, 0]
        if np.dot(vector, formal_pc1) < 0:
            vector *= -1
        boot_loadings[rep, :] = vector
    pca_boot_rows = []
    loading_summary_rows = []
    major = np.abs(formal_pc1) >= 0.35
    stable_major_count = 0
    for ai, axis in enumerate(AXES):
        same_sign = boot_loadings[:, ai] * formal_pc1[ai] > 0
        sign_stability = float(np.mean(same_sign))
        if major[ai] and sign_stability >= 0.80:
            stable_major_count += 1
        loading_summary_rows.append({
            "axis": axis, "loading": formal_pc1[ai], "bootstrap_median": np.median(boot_loadings[:, ai]), "bootstrap_ci_low": np.quantile(boot_loadings[:, ai], 0.025), "bootstrap_ci_high": np.quantile(boot_loadings[:, ai], 0.975),
            "formal_sign": "POSITIVE" if formal_pc1[ai] > 0 else "NEGATIVE", "sign_stability_proportion": sign_stability, "stable_loading_ge_0_80": sign_stability >= 0.80,
            "major_loading_abs_ge_0_35": bool(major[ai]), "sign_anchor_axis": anchor_axis, "orientation_rule": orientation_rule,
        })
        for rep in range(B):
            pca_boot_rows.append({"replicate": rep + 1, "axis": axis, "loading": boot_loadings[rep, ai], "aligned_to_formal_PC1": True})
    loading_summary = pd.DataFrame(loading_summary_rows)
    write_tsv(loading_summary, "03_realization_axis/CGGA325_PC1_LOADINGS.tsv")
    write_tsv(pd.DataFrame(pca_boot_rows), "03_realization_axis/CGGA325_PC1_BOOTSTRAP.tsv")
    patient_scores = discovery[["cohort", "patient_id", "bulk_mes_z"]].copy()
    patient_scores["REALIZATION_PC1_DISCOVERY"] = scores_all[:, 0]
    write_tsv(patient_scores, "03_realization_axis/CGGA325_PATIENT_REALIZATION_SCORES.tsv")

    gate1 = variance_fraction[0] >= 0.35
    gate2 = int(major.sum()) >= 2
    gate3 = stable_major_count >= 2
    discovery_pass = bool(gate1 and gate2 and gate3)
    discovery_gate = pd.DataFrame([{
        "PC1_variance_fraction": variance_fraction[0], "gate1_variance_ge_0_35": gate1, "n_abs_loading_ge_0_35": int(major.sum()), "gate2_at_least_two_major": gate2,
        "n_major_with_sign_stability_ge_0_80": stable_major_count, "gate3_at_least_two_stable_major": gate3,
        "discovery_gate_status": "DISCOVERY_REALIZATION_AXIS_PASS" if discovery_pass else "NO_STABLE_DOMINANT_REALIZATION_AXIS_IN_DISCOVERY",
        "sign_anchor_axis": anchor_axis, "orientation_rule": orientation_rule,
    }])
    write_tsv(discovery_gate, "03_realization_axis/CGGA325_DISCOVERY_GATE.tsv")

    locked = pd.DataFrame({"axis_order": np.arange(1, 6), "axis": AXES, "loading": formal_pc1, "sign_anchor_axis": anchor_axis, "orientation_rule": orientation_rule, "discovery_gate_status": discovery_gate.discovery_gate_status.iloc[0]})
    locked_path = RUN / "03_realization_axis/CGGA325_PC1_LOADINGS_LOCKED.tsv"
    locked.to_csv(locked_path, sep="\t", index=False, lineterminator="\n")
    locked_sha = sha256(locked_path)
    (RUN / "03_realization_axis/CGGA325_PC1_LOADINGS_LOCKED.tsv.sha256").write_text(f"{locked_sha}  {locked_path.name}\n", encoding="utf-8")

    completed = datetime.now(timezone.utc)
    log = {
        "started_at_utc": started.isoformat().replace("+00:00", "Z"),
        "completed_at_utc": completed.isoformat().replace("+00:00", "Z"),
        "protocol_sha256": protocol_sha,
        "freeze_validation": "PASS",
        "v1_6_checksum_files_verified": checked,
        "correlation_bootstraps_per_cohort": B,
        "global_permutations_per_cohort": B,
        "pca_bootstraps": B,
        "level1_gate": "PASS" if level1 else "FAIL",
        "discovery_gate": discovery_gate.discovery_gate_status.iloc[0],
        "locked_loading_sha256": locked_sha,
    }
    (RUN / "logs/residual_structure_discovery.log").write_text(json.dumps(log, indent=2) + "\n", encoding="utf-8")
    session = {"python": platform.python_version(), "platform": platform.platform(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "script_sha256": sha256(Path(__file__)), "biweight_status": "NOT_RUN_DEPENDENCY_NOT_PRECERTIFIED"}
    (RUN / "session_info/RESIDUAL_STRUCTURE_DISCOVERY_SESSION_INFO.json").write_text(json.dumps(session, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(log, indent=2))


if __name__ == "__main__":
    main()
