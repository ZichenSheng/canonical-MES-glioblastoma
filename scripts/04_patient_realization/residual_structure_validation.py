#!/usr/bin/env python3
"""Formal locked CGGA693 projection, transfer gate, and orthogonality QA."""

from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
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


def pca(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    centered = x - x.mean(axis=0, keepdims=True)
    _, singular, vt = np.linalg.svd(centered, full_matrices=False)
    eigenvalues = singular**2 / (x.shape[0] - 1)
    loadings = vt.T
    scores = centered @ loadings
    return eigenvalues, loadings, scores


def main() -> None:
    started = datetime.now(timezone.utc)
    seeds = json.loads((RUN / "00_governance/V1_7_RANDOM_SEEDS.json").read_text())["seeds"]
    locked_path = RUN / "03_realization_axis/CGGA325_PC1_LOADINGS_LOCKED.tsv"
    sidecar = RUN / "03_realization_axis/CGGA325_PC1_LOADINGS_LOCKED.tsv.sha256"
    expected_sha = sidecar.read_text().split()[0]
    observed_sha = sha256(locked_path)
    if expected_sha != observed_sha:
        raise RuntimeError("Locked discovery loading SHA256 mismatch")
    locked = pd.read_csv(locked_path, sep="\t").sort_values("axis_order")
    if list(locked.axis) != AXES:
        raise RuntimeError("Locked axis order mismatch")
    discovery_loading = locked.loading.to_numpy(float)
    discovery_gate_status = locked.discovery_gate_status.iloc[0]
    discovery_pass = discovery_gate_status == "DISCOVERY_REALIZATION_AXIS_PASS"
    anchor_axis = locked.sign_anchor_axis.iloc[0]
    anchor_index = AXES.index(anchor_axis)

    residual = pd.read_csv(V16 / "02_fingerprint/MES_CONDITIONED_FINGERPRINT.tsv", sep="\t")
    master = pd.read_csv(V16 / "01_reference/PATIENT_REALIZATION_MASTER.tsv", sep="\t")
    scores325 = pd.read_csv(RUN / "03_realization_axis/CGGA325_PATIENT_REALIZATION_SCORES.tsv", sep="\t")
    data693 = residual.loc[residual.cohort == "CGGA693"].copy().sort_values("patient_id").reset_index(drop=True)
    if len(data693) != 108 or data693[COLS].isna().any().any():
        raise RuntimeError("CGGA693 validation matrix mismatch")
    x693 = data693[COLS].to_numpy(float)
    projected = x693 @ discovery_loading
    projected_table = data693[["cohort", "patient_id", "bulk_mes_z"]].copy()
    projected_table["REALIZATION_PC1_PROJECTED_693"] = projected
    projected_table["locked_loading_sha256"] = observed_sha
    projected_table["axis_status"] = "LOCKED_DIAGNOSTIC_PROJECTION" if not discovery_pass else "LOCKED_VALIDATION_PROJECTION"
    write_tsv(projected_table, "04_cross_cohort_validation/CGGA693_PROJECTED_REALIZATION.tsv")

    eigen, val_loadings, val_scores = pca(x693)
    val_pc1 = val_loadings[:, 0].copy()
    # The discovery anchor axis fixes validation sign without narrative intervention.
    if val_pc1[anchor_index] < 0:
        val_pc1 *= -1
        val_loadings[:, 0] *= -1
        val_scores[:, 0] *= -1
    variance = eigen / eigen.sum()
    independent = pd.DataFrame({
        "component": [f"PC{i}" for i in range(1, 6)],
        "eigenvalue": eigen,
        "variance_fraction": variance,
        "cumulative_variance_fraction": np.cumsum(variance),
        "validation_role": ["INDEPENDENT_PC1_DIAGNOSTIC"] + ["DESCRIPTIVE_QA_ONLY"] * 4,
    })
    for ai, axis in enumerate(AXES):
        independent[f"loading_{axis}"] = val_loadings[ai, :]
    write_tsv(independent, "04_cross_cohort_validation/CGGA693_INDEPENDENT_PCA.tsv")
    val_loading_table = pd.DataFrame({"axis": AXES, "discovery_loading_locked": discovery_loading, "validation_independent_PC1_loading": val_pc1})
    write_tsv(val_loading_table, "04_cross_cohort_validation/CGGA693_INDEPENDENT_PC1_LOADINGS.tsv")

    cosine = float(np.dot(discovery_loading, val_pc1) / (np.linalg.norm(discovery_loading) * np.linalg.norm(val_pc1)))
    loading_r, loading_p = pearsonr(discovery_loading, val_pc1)
    total_validation_variance = float(np.var(x693, axis=0, ddof=1).sum())
    projected_variance = float(np.var(projected, ddof=1))
    projected_fraction = projected_variance / total_validation_variance

    relation_rows = []
    for ai, axis in enumerate(AXES):
        sr, sp = spearmanr(projected, x693[:, ai])
        pr, pp = pearsonr(projected, x693[:, ai])
        relation_rows.append({
            "cohort": "CGGA693", "axis": axis, "discovery_loading": discovery_loading[ai],
            "projected_score_spearman": sr, "spearman_p": sp, "projected_score_pearson": pr, "pearson_p": pp,
            "direction_matches_discovery_loading": bool(np.sign(sr) == np.sign(discovery_loading[ai])),
            "absolute_discovery_loading_rank": int(pd.Series(np.abs(discovery_loading)).rank(method="first", ascending=False).iloc[ai]),
        })
    relations = pd.DataFrame(relation_rows)
    write_tsv(relations, "04_cross_cohort_validation/PROJECTED_AXIS_COMPONENT_RELATIONS.tsv")

    top_two_indices = np.argsort(-np.abs(discovery_loading))[:2]
    top_two_axes = [AXES[i] for i in top_two_indices]
    top_two_concordant = bool(relations.loc[relations.axis.isin(top_two_axes), "direction_matches_discovery_loading"].all())
    gate_cosine = cosine >= 0.70
    gate_variance = projected_fraction >= 0.25
    gate_top_two = top_two_concordant
    transfer_pass = bool(discovery_pass and gate_cosine and gate_variance and gate_top_two)
    concordance = pd.DataFrame([{
        "locked_loading_sha256": observed_sha, "discovery_sign_anchor_axis": anchor_axis,
        "cosine_similarity": cosine, "pearson_correlation_of_five_loadings": loading_r, "pearson_p_descriptive_n5": loading_p,
        "validation_independent_PC1_variance_fraction": variance[0], "projected_variance": projected_variance,
        "total_validation_residual_variance": total_validation_variance, "projected_variance_fraction": projected_fraction,
        "interpretation_boundary": "FIVE_LOADING_EFFECT_SIZES; VALIDATION_PCA_DID_NOT_MODIFY_LOCKED_AXIS",
    }])
    write_tsv(concordance, "04_cross_cohort_validation/LOADING_CONCORDANCE.tsv")
    gate = pd.DataFrame([{
        "discovery_gate_status": discovery_gate_status, "gate1_discovery_pass": discovery_pass,
        "cosine_similarity": cosine, "gate2_cosine_ge_0_70": gate_cosine,
        "projected_variance_fraction": projected_fraction, "gate3_projected_fraction_ge_0_25": gate_variance,
        "top_two_discovery_axes": ";".join(top_two_axes), "top_two_directional_concordance": top_two_concordant, "gate4_top_two_direction_match": gate_top_two,
        "cross_cohort_gate_status": "CROSS_COHORT_REALIZATION_AXIS_SUPPORTED" if transfer_pass else "REALIZATION_AXIS_NOT_STABLY_TRANSFERABLE",
        "interpretation_status": "VALIDATED_AXIS" if transfer_pass else "DIAGNOSTIC_PROJECTION_ONLY",
    }])
    write_tsv(gate, "04_cross_cohort_validation/CROSS_COHORT_REALIZATION_GATE.tsv")

    # Module J: patient-bootstrap orthogonality QA for discovery and locked validation scores.
    orth_rows = []
    for cohort_index, cohort in enumerate(["CGGA325", "CGGA693"]):
        if cohort == "CGGA325":
            table = scores325[["patient_id", "bulk_mes_z", "REALIZATION_PC1_DISCOVERY"]].copy()
            score = table.REALIZATION_PC1_DISCOVERY.to_numpy(float)
            score_name = "REALIZATION_PC1_DISCOVERY"
        else:
            table = projected_table[["patient_id", "bulk_mes_z", "REALIZATION_PC1_PROJECTED_693"]].copy()
            score = table.REALIZATION_PC1_PROJECTED_693.to_numpy(float)
            score_name = "REALIZATION_PC1_PROJECTED_693"
        mes = table.bulk_mes_z.to_numpy(float)
        rng = np.random.default_rng(seeds["orthogonality_bootstrap"] + cohort_index * 10000)
        boot_s = np.empty(B, float); boot_p = np.empty(B, float)
        for rep in range(B):
            idx = rng.integers(0, len(mes), len(mes))
            boot_s[rep] = spearmanr(score[idx], mes[idx]).statistic
            boot_p[rep] = pearsonr(score[idx], mes[idx]).statistic
        for method, observed, pvalue, boot_values in [
            ("SPEARMAN", spearmanr(score, mes).statistic, spearmanr(score, mes).pvalue, boot_s),
            ("PEARSON", pearsonr(score, mes).statistic, pearsonr(score, mes).pvalue, boot_p),
        ]:
            orth_rows.append({
                "cohort": cohort, "score": score_name, "method": method, "estimate": observed,
                "bootstrap_ci_low": np.quantile(boot_values, 0.025), "bootstrap_ci_high": np.quantile(boot_values, 0.975),
                "nominal_p": pvalue, "bootstrap_replicates": B,
                "orthogonality_status": "RESIDUAL_AXIS_RETAINS_MATERIAL_MES_ASSOCIATION" if method == "SPEARMAN" and abs(observed) > 0.20 else "NO_MATERIAL_MES_ASSOCIATION",
            })
    write_tsv(pd.DataFrame(orth_rows), "04_cross_cohort_validation/ORTHOGONALITY_QA.tsv")

    # Module K is strictly conditional on the complete cross-cohort gate.
    if transfer_pass:
        pair_data = pd.read_csv(V16 / "03_iso_score/ALL_PATIENT_PAIR_DISTANCES.tsv", sep="\t")
        hero = pd.read_csv(V16 / "03_iso_score/ISO_MES_HERO_PAIRS.tsv", sep="\t")
        score_maps = {
            "CGGA325": dict(zip(scores325.patient_id, scores325.REALIZATION_PC1_DISCOVERY)),
            "CGGA693": dict(zip(projected_table.patient_id, projected_table.REALIZATION_PC1_PROJECTED_693)),
        }
        link_rows = []
        for cohort in ["CGGA325", "CGGA693"]:
            pairs = pair_data.loc[(pair_data.cohort == cohort) & (pair_data.delta_mes <= 0.10)].copy()
            smap = score_maps[cohort]
            pairs["abs_delta_realization"] = (pairs.patient_i.map(smap) - pairs.patient_j.map(smap)).abs()
            sr = spearmanr(pairs.raw_fingerprint_distance, pairs.abs_delta_realization)
            pr = pearsonr(pairs.raw_fingerprint_distance, pairs.abs_delta_realization)
            h = hero.loc[hero.cohort == cohort].iloc[0]
            hero_delta = abs(smap[h.patient_i] - smap[h.patient_j])
            link_rows.append({
                "cohort": cohort, "status": "EVALUATED_AFTER_TRANSFER_GATE", "n_iso_pairs": len(pairs),
                "median_abs_delta_realization": pairs.abs_delta_realization.median(), "P90_abs_delta_realization": pairs.abs_delta_realization.quantile(0.9),
                "hero_patient_i": h.patient_i, "hero_patient_j": h.patient_j, "hero_abs_delta_realization": hero_delta,
                "distance_relation_spearman": sr.statistic, "spearman_p_descriptive": sr.pvalue,
                "distance_relation_pearson": pr.statistic, "pearson_p_descriptive": pr.pvalue,
                "interpretation_boundary": "PAIRWISE_DESCRIPTIVE_LINK_NOT_PERCENT_DIVERGENCE_EXPLAINED",
            })
        link = pd.DataFrame(link_rows)
    else:
        link = pd.DataFrame([{
            "cohort": "CGGA325;CGGA693", "status": "NOT_RUN_CROSS_COHORT_REALIZATION_GATE_FAILED",
            "n_iso_pairs": np.nan, "median_abs_delta_realization": np.nan, "P90_abs_delta_realization": np.nan,
            "hero_patient_i": np.nan, "hero_patient_j": np.nan, "hero_abs_delta_realization": np.nan,
            "distance_relation_spearman": np.nan, "spearman_p_descriptive": np.nan,
            "distance_relation_pearson": np.nan, "pearson_p_descriptive": np.nan,
            "interpretation_boundary": "LOCKED_PC1_REMAINS_DIAGNOSTIC_AND_IS_NOT_USED_TO_INTERPRET_V1_6_PAIR_DIVERGENCE",
        }])
    write_tsv(link, "04_cross_cohort_validation/ISO_MES_REALIZATION_LINK.tsv")

    completed = datetime.now(timezone.utc)
    log = {
        "started_at_utc": started.isoformat().replace("+00:00", "Z"),
        "completed_at_utc": completed.isoformat().replace("+00:00", "Z"),
        "locked_loading_sha256_verified": observed_sha,
        "discovery_gate": discovery_gate_status,
        "cosine_similarity": cosine,
        "projected_variance_fraction": projected_fraction,
        "top_two_axes": top_two_axes,
        "top_two_directional_concordance": top_two_concordant,
        "cross_cohort_gate": gate.cross_cohort_gate_status.iloc[0],
        "iso_mes_link_status": link.status.iloc[0],
    }
    (RUN / "logs/locked_validation.log").write_text(json.dumps(log, indent=2) + "\n", encoding="utf-8")
    session = {"python": platform.python_version(), "platform": platform.platform(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "script_sha256": sha256(Path(__file__))}
    (RUN / "session_info/LOCKED_VALIDATION_SESSION_INFO.json").write_text(json.dumps(session, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(log, indent=2))


if __name__ == "__main__":
    main()
