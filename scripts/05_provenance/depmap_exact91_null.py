#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import re
from collections import defaultdict
import os
from pathlib import Path

DATA_ROOT = Path(os.environ.get("GBM_MES_DATA_ROOT", "data"))
RESULT_ROOT = Path(os.environ.get("GBM_MES_OUTPUT_ROOT", "results"))
REPO_ROOT = Path(os.environ.get("GBM_MES_REPO_ROOT", Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd
from scipy.spatial.distance import cosine
from scipy.stats import rankdata, spearmanr

ROOT = RESULT_ROOT / "depmap_null"
V21 = DATA_ROOT / "prepared" / "depmap_controls"
V2 = RESULT_ROOT / "depmap"
DATA = DATA_ROOT / "depmap"
SEED = 20260812
NBOOT = 2000
NPERM = 5000
NRANDOM = 5000
MIN_PAIRS = 30
GENE_RE = re.compile(r"\s+\(\d+\)$")


def clean_gene(x: str) -> str:
    return GENE_RE.sub("", str(x)).upper()


def greedy_caliper(a: np.ndarray, b: np.ndarray, caliper: float) -> list[tuple[int, int]]:
    candidates = sorted((abs(a[i] - b[j]), i, j) for i in range(len(a)) for j in range(len(b)) if abs(a[i] - b[j]) <= caliper)
    used_a: set[int] = set()
    used_b: set[int] = set()
    out: list[tuple[int, int]] = []
    for _, i, j in candidates:
        if i not in used_a and j not in used_b:
            out.append((i, j)); used_a.add(i); used_b.add(j)
    return out


def match_ids(ng_ids: list[str], tr_ids: list[str], score: pd.Series) -> tuple[list[tuple[str, str]], float]:
    a = score.loc[ng_ids].to_numpy(float)
    b = score.loc[tr_ids].to_numpy(float)
    caliper = float(.2 * np.std(np.r_[a, b], ddof=1))
    return [(ng_ids[i], tr_ids[j]) for i, j in greedy_caliper(a, b, caliper)], caliper


def corr(x: np.ndarray, y: np.ndarray) -> float:
    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(spearmanr(x, y).statistic)


def pair_metrics(pairs: list[tuple[str, str]], z: pd.DataFrame, genes: list[str], score: pd.Series) -> pd.DataFrame:
    rows = []
    p = len(genes)
    for i, (ng, tr) in enumerate(pairs, 1):
        x = z.loc[ng, genes].to_numpy(float)
        y = z.loc[tr, genes].to_numpy(float)
        rd = float(np.mean(np.abs(rankdata(x) - rankdata(y))) / (p - 1))
        rows.append({
            "pair_id": i, "NextGen_ModelID": ng, "traditional_ModelID": tr,
            "NextGen_score": score.loc[ng], "traditional_score": score.loc[tr],
            "absolute_score_difference": abs(score.loc[ng] - score.loc[tr]),
            "contribution_correlation": corr(x, y),
            "cosine_similarity": float(1 - cosine(x, y)) if np.linalg.norm(x) and np.linalg.norm(y) else np.nan,
            "rank_similarity": 1 - rd,
            "contribution_distance": float(np.sqrt(np.mean(((x - y) / p) ** 2))),
        })
    return pd.DataFrame(rows)


def bootstrap_mean(x: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    reps = np.empty(NBOOT)
    for i in range(NBOOT):
        reps[i] = np.mean(rng.choice(x, len(x), replace=True))
    return float(np.quantile(reps, .025)), float(np.quantile(reps, .975))


def label_permutation(ids: list[str], n_ng: int, score: pd.Series, z: pd.DataFrame, genes: list[str], rng: np.random.Generator) -> np.ndarray:
    ids = np.asarray(ids)
    out = np.empty(NPERM)
    for i in range(NPERM):
        sh = rng.permutation(ids)
        aa, bb = sh[:n_ng].tolist(), sh[n_ng:].tolist()
        pairs, _ = match_ids(aa, bb, score)
        out[i] = pair_metrics(pairs, z, genes, score).contribution_correlation.mean()
    return out


def qbins(values: pd.Series, q: int) -> pd.Series:
    return pd.qcut(values.rank(method="first"), q=q, labels=False, duplicates="drop").astype(int)


def empirical_lower(null: np.ndarray, observed: float) -> float:
    null = null[np.isfinite(null)]
    return float((1 + np.sum(null <= observed)) / (len(null) + 1))


def coex_matrix(x: np.ndarray) -> np.ndarray:
    ranks = rankdata(x, axis=0, method="average")
    ranks -= ranks.mean(axis=0)
    den = np.sqrt(np.sum(ranks * ranks, axis=0))
    ranks /= den
    return ranks.astype(np.float32)


def coex_summary_from_norm(norm: np.ndarray, indices: list[int]) -> dict[str, float]:
    x = norm[:, indices]
    c = x.T @ x
    vals = np.abs(c[np.triu_indices(len(indices), 1)]).astype(float)
    q25, q50, q75 = np.quantile(vals, [.25, .50, .75])
    return {
        "mean_abs_rho": float(np.mean(vals)), "median_abs_rho": float(q50),
        "q25_abs_rho": float(q25), "q75_abs_rho": float(q75), "IQR_abs_rho": float(q75 - q25),
        "fraction_abs_rho_gt_0.2": float(np.mean(vals > .2)), "fraction_abs_rho_gt_0.4": float(np.mean(vals > .4)),
    }


def candidate_pool(features: pd.DataFrame, mes: set[str]) -> dict[tuple[int, int, int], list[str]]:
    pools: dict[tuple[int, int, int], list[str]] = defaultdict(list)
    for gene, row in features.loc[~features.index.isin(mes)].iterrows():
        pools[(int(row.mean_bin), int(row.var_bin), int(row.detect_bin))].append(gene)
    return pools


def available_from_keys(pools: dict, keys: list[tuple[int, int, int]], used: set[str]) -> list[str]:
    out = []
    for key in keys:
        out.extend(g for g in pools.get(key, []) if g not in used)
    return list(dict.fromkeys(out))


def generate_signature(template_genes: list[str], features: pd.DataFrame, pools: dict, rng: np.random.Generator) -> tuple[list[str], list[dict], bool]:
    used: set[str] = set()
    selected: list[str] = []
    check: list[dict] = []
    for template in template_genes:
        row = features.loc[template]
        m, v, d = int(row.mean_bin), int(row.var_bin), int(row.detect_bin)
        levels = [
            ("EXACT_BIN", [(m, v, d)]),
            ("L1_ADJACENT_DETECTION", [(m, v, dd) for dd in (d - 1, d + 1) if 0 <= dd <= 4]),
            ("L2_ADJACENT_VARIANCE", [(m, vv, d) for vv in (v - 1, v + 1) if 0 <= vv <= 9]),
            ("L3_ADJACENT_MEAN_VARIANCE", [(mm, vv, d) for mm in (m - 1, m + 1) for vv in (v - 1, v + 1) if 0 <= mm <= 9 and 0 <= vv <= 9]),
        ]
        pool, level = [], "FAILED"
        for name, keys in levels:
            pool = available_from_keys(pools, keys, used)
            if pool:
                level = name; break
        if not pool:
            return selected, check, False
        gene = str(rng.choice(pool))
        used.add(gene); selected.append(gene)
        check.append({
            "position": len(selected), "template_MES_gene": template, "random_gene": gene, "match_level": level,
            "mean_expression_distance": abs(features.loc[gene, "mean"] - row["mean"]),
            "variance_distance": abs(features.loc[gene, "variance"] - row["variance"]),
            "detection_distance": abs(features.loc[gene, "detection"] - row["detection"]),
        })
    return selected, check, len(selected) == 91 and len(set(selected)) == 91


def main() -> None:
    rng = np.random.default_rng(SEED)
    sig = pd.read_csv(V2 / "02_signature_reference/signature_registry.tsv", sep="\t")
    scores_v2 = pd.read_csv(V2 / "03_expression_scores/model_level_scores.tsv", sep="\t").set_index("ModelID")
    ng_ids = scores_v2.index[scores_v2.NextGen_or_traditional.eq("NextGen")].tolist()
    tr_ids = scores_v2.index[scores_v2.NextGen_or_traditional.eq("traditional")].tolist()
    assert (len(ng_ids), len(tr_ids)) == (70, 55)

    headers = pd.read_csv(DATA / "next_gen_expression.csv", nrows=0).columns.tolist()
    dtypes = {c: np.float32 for c in headers[1:]}; dtypes[headers[0]] = str
    ng_raw = pd.read_csv(DATA / "next_gen_expression.csv", index_col=0, dtype=dtypes)
    tr_raw = pd.read_csv(DATA / "traditional_expression.csv", index_col=0, dtype=dtypes)
    raw = pd.concat([ng_raw, tr_raw]); raw.index = raw.index.astype(str)
    symbol_to_col: dict[str, str] = {}
    for col in raw.columns:
        symbol_to_col.setdefault(clean_gene(col), col)
    cols = list(symbol_to_col.values())
    unique_raw = raw[cols].copy(); unique_raw.columns = [clean_gene(c) for c in cols]
    mu = unique_raw.mean(); sd = unique_raw.std(ddof=1).replace(0, np.nan)
    strict_ids = ng_ids + tr_ids
    z = ((unique_raw.loc[strict_ids] - mu) / sd).astype(np.float32)

    reference = sig.loc[sig.signature_id.eq("CANONICAL_MES95"), "gene"].astype(str).str.upper().drop_duplicates().tolist()
    measurable = [g for g in reference if g in unique_raw.columns]
    missing = [g for g in reference if g not in unique_raw.columns]
    assert len(reference) == 95 and len(measurable) == 91 and len(missing) == 4
    auth_rows = []
    for gene in reference:
        measured = gene in unique_raw.columns
        auth_rows.append({"gene": gene, "reference_status": "CANONICAL_MES95_V1_4_FROZEN", "measured_nextgen": measured, "measured_traditional": measured, "included_primary": measured, "reason_if_missing": "NOT_PRESENT_IN_CERTIFIED_EXPRESSION_COLUMNS_NO_IMPUTATION" if not measured else ""})
    auth_df = pd.DataFrame(auth_rows)
    auth_df.to_csv(ROOT / "01_measurable_mes95/mes95_reference_95.tsv", sep="\t", index=False)
    auth_df.loc[auth_df.included_primary].to_csv(ROOT / "01_measurable_mes95/mes95_measurable_91.tsv", sep="\t", index=False)
    auth_df.loc[~auth_df.included_primary].to_csv(ROOT / "01_measurable_mes95/mes95_missing_4.tsv", sep="\t", index=False)

    # Canonical 91-gene reproduction.
    mes_score = z[measurable].mean(axis=1)
    reconstruction_error = float(np.max(np.abs(mes_score - scores_v2["score_z__CANONICAL_MES95"])))
    assert reconstruction_error < 1e-5
    mes_pairs, mes_caliper = match_ids(ng_ids, tr_ids, mes_score)
    mes_pair_metrics = pair_metrics(mes_pairs, z, measurable, mes_score)
    mes_corr = float(mes_pair_metrics.contribution_correlation.mean())
    mes_ci = bootstrap_mean(mes_pair_metrics.contribution_correlation.to_numpy(float), rng)
    mes_perm = label_permutation(strict_ids, len(ng_ids), mes_score, z, measurable, rng)
    mes_perm_p = empirical_lower(mes_perm, mes_corr)
    absdiff = mes_pair_metrics.absolute_score_difference.to_numpy(float)
    pa = np.array([mes_score.loc[x] for x, _ in mes_pairs]); pb = np.array([mes_score.loc[y] for _, y in mes_pairs])
    post_smd = float((pa.mean() - pb.mean()) / math.sqrt((pa.var(ddof=1) + pb.var(ddof=1)) / 2))
    observed = pd.DataFrame([{
        "reference_gene_n": 95, "measurable_gene_n": 91, "matched_pair_n": len(mes_pairs), "caliper": mes_caliper,
        "absolute_score_difference_mean": absdiff.mean(), "absolute_score_difference_median": np.median(absdiff),
        "absolute_score_difference_IQR": np.quantile(absdiff, .75) - np.quantile(absdiff, .25), "absolute_score_difference_max": absdiff.max(),
        "post_match_SMD": post_smd, "contribution_correlation": mes_corr,
        "bootstrap_CI_low": mes_ci[0], "bootstrap_CI_high": mes_ci[1], "label_permutation_P": mes_perm_p,
        "score_reconstruction_max_abs_error": reconstruction_error,
    }])
    observed.to_csv(ROOT / "02_exact91_random_null/mes95_exact91_observed.tsv", sep="\t", index=False)
    mes_pair_metrics.to_csv(ROOT / "02_exact91_random_null/mes95_exact91_pairs.tsv", sep="\t", index=False)
    pd.DataFrame({"iteration": np.arange(1, NPERM + 1), "mean_contribution_correlation": mes_perm}).to_csv(ROOT / "02_exact91_random_null/mes95_label_permutation_null.tsv.gz", sep="\t", index=False, compression="gzip")

    # Feature-matched exact 91-gene candidates, generated without outcome access.
    features = pd.DataFrame({"mean": mu, "variance": unique_raw.var(ddof=1), "detection": (unique_raw > 0).mean()}).replace([np.inf, -np.inf], np.nan).dropna()
    features["mean_bin"] = qbins(features["mean"], 10)
    features["var_bin"] = qbins(features["variance"], 10)
    features["detect_bin"] = qbins(features["detection"], 5)
    pools = candidate_pool(features, set(reference))
    random_sets: list[list[str]] = []
    membership_rows: list[dict] = []
    generation_rows = []
    for sid in range(1, NRANDOM + 1):
        genes, check, success = generate_signature(measurable, features, pools, rng)
        random_sets.append(genes)
        for row in check:
            membership_rows.append({"signature_id": sid, **row})
        generation_rows.append({
            "signature_id": sid, "n_genes": len(genes), "matching_success": success,
            "n_exact_bin": sum(r["match_level"] == "EXACT_BIN" for r in check),
            "n_near_bin": sum(r["match_level"] != "EXACT_BIN" for r in check),
            "mean_expression_distance": np.mean([r["mean_expression_distance"] for r in check]) if check else np.nan,
            "variance_distance": np.mean([r["variance_distance"] for r in check]) if check else np.nan,
            "detection_distance": np.mean([r["detection_distance"] for r in check]) if check else np.nan,
        })
    membership = pd.DataFrame(membership_rows)
    membership.to_csv(ROOT / "02_exact91_random_null/exact91_random_membership.tsv.gz", sep="\t", index=False, compression="gzip")

    # Coexpression matrices are outcome-blind: standardize each gene within context first,
    # concatenate residualized contexts, then rank genes across the combined model universe.
    ngz = z.loc[ng_ids]
    trz = z.loc[tr_ids]
    residual = pd.concat([(ngz - ngz.mean()) / ngz.std(ddof=1), (trz - trz.mean()) / trz.std(ddof=1)])
    combined_norm = coex_matrix(residual.to_numpy(float))
    ng_norm = coex_matrix(ngz.to_numpy(float))
    tr_norm = coex_matrix(trz.to_numpy(float))
    gene_index = {g: i for i, g in enumerate(z.columns)}
    mes_idx = [gene_index[g] for g in measurable]
    mes_combined = coex_summary_from_norm(combined_norm, mes_idx)
    mes_ng = coex_summary_from_norm(ng_norm, mes_idx)
    mes_tr = coex_summary_from_norm(tr_norm, mes_idx)
    pd.DataFrame([{"context": "combined_within_context_standardized", "model_n": len(strict_ids), **mes_combined}, {"context": "NextGen", "model_n": len(ng_ids), **mes_ng}, {"context": "traditional", "model_n": len(tr_ids), **mes_tr}]).to_csv(ROOT / "03_coexpression_null/mes95_coexpression_summary.tsv", sep="\t", index=False)

    random_rows = []
    for sid, genes in enumerate(random_sets, 1):
        base = generation_rows[sid - 1]
        if not base["matching_success"]:
            random_rows.append({**base, "n_matched_pairs": 0, "eligible_signature": False}); continue
        score = z[genes].mean(axis=1)
        pairs, caliper = match_ids(ng_ids, tr_ids, score)
        pm = pair_metrics(pairs, z, genes, score)
        cx = coex_summary_from_norm(combined_norm, [gene_index[g] for g in genes])
        random_rows.append({
            **base, "caliper": caliper, "n_matched_pairs": len(pairs), "eligible_signature": len(pairs) >= MIN_PAIRS,
            "contribution_correlation": pm.contribution_correlation.mean(), "cosine_similarity": pm.cosine_similarity.mean(),
            "rank_similarity": pm.rank_similarity.mean(), "contribution_distance": pm.contribution_distance.mean(),
            **{f"coexpression_{k}": v for k, v in cx.items()},
        })
    random_df = pd.DataFrame(random_rows)
    random_df.to_csv(ROOT / "02_exact91_random_null/exact91_random_null.tsv", sep="\t", index=False)
    eligible = random_df.loc[random_df.matching_success.eq(True) & random_df.eligible_signature.eq(True)].copy()
    null_corr = eligible.contribution_correlation.to_numpy(float)
    exact_p = empirical_lower(null_corr, mes_corr)
    divergence_percentile = float(100 * np.mean(null_corr >= mes_corr))
    z_position = float((mes_corr - null_corr.mean()) / null_corr.std(ddof=1))
    exact_decision = "EXACT_GENE_COUNT_NULL_PASS" if divergence_percentile >= 95 and exact_p < .05 else "EXACT_GENE_COUNT_NULL_FAIL"
    exact_summary = pd.DataFrame([{
        "random_signatures_requested": NRANDOM, "matching_success_n": int(random_df.matching_success.sum()),
        "eligible_signatures_n": len(eligible), "minimum_pairs": MIN_PAIRS, "gene_n": 91,
        "MES95_matched_pair_n": len(mes_pairs), "MES95_contribution_correlation": mes_corr,
        "random_null_mean": null_corr.mean(), "random_null_SD": null_corr.std(ddof=1),
        "MES95_divergence_percentile": divergence_percentile, "empirical_P": exact_p,
        "z_position_similarity_scale": z_position, "effect_MES95_minus_null_mean": mes_corr - null_corr.mean(),
        "decision": exact_decision,
    }])
    exact_summary.to_csv(ROOT / "02_exact91_random_null/exact91_random_summary.tsv", sep="\t", index=False)

    # Coexpression-aware selection from outcome-blind exact91 candidates.
    selection_rows = []
    target_med, target_iqr = mes_combined["median_abs_rho"], mes_combined["IQR_abs_rho"]
    for tol in (.10, .15, .20):
        mask = eligible.coexpression_median_abs_rho.between(target_med * (1 - tol), target_med * (1 + tol)) & eligible.coexpression_IQR_abs_rho.between(target_iqr * (1 - tol), target_iqr * (1 + tol))
        selection_rows.append({"relative_tolerance": tol, "eligible_signature_n": int(mask.sum()), "target_median_abs_rho": target_med, "target_IQR_abs_rho": target_iqr})
    counts = {r["relative_tolerance"]: r["eligible_signature_n"] for r in selection_rows}
    chosen = next((t for t in (.10, .15, .20) if counts[t] >= 1000), None)
    if chosen is None:
        chosen = next((t for t in (.10, .15, .20) if counts[t] >= 500), .20)
    coex_mask = eligible.coexpression_median_abs_rho.between(target_med * (1 - chosen), target_med * (1 + chosen)) & eligible.coexpression_IQR_abs_rho.between(target_iqr * (1 - chosen), target_iqr * (1 + chosen))
    coex_null = eligible.loc[coex_mask].copy()
    effective_n = len(coex_null)
    coex_corr = coex_null.contribution_correlation.to_numpy(float)
    coex_p = empirical_lower(coex_corr, mes_corr) if effective_n else np.nan
    coex_percentile = float(100 * np.mean(coex_corr >= mes_corr)) if effective_n else np.nan
    if effective_n < 500:
        coex_decision = "COEXPRESSION_NULL_LOW_POWER"
    elif coex_percentile >= 95 and coex_p < .05:
        coex_decision = "COEXPRESSION_AWARE_NULL_PASS"
    else:
        coex_decision = "COEXPRESSION_AWARE_NULL_FAIL"
    coex_null.to_csv(ROOT / "03_coexpression_null/coexpression_null.tsv", sep="\t", index=False)
    coex_summary = pd.DataFrame([{
        "candidate_signatures_n": len(eligible), "effective_random_signatures_n": effective_n,
        "chosen_relative_tolerance": chosen, "MES95_median_abs_rho": target_med, "MES95_IQR_abs_rho": target_iqr,
        "MES95_contribution_correlation": mes_corr, "coexpression_null_mean": np.mean(coex_corr) if effective_n else np.nan,
        "coexpression_null_SD": np.std(coex_corr, ddof=1) if effective_n > 1 else np.nan,
        "MES95_divergence_percentile": coex_percentile, "empirical_P": coex_p, "decision": coex_decision,
        "selection_did_not_use_outcome": True,
    }])
    coex_summary.to_csv(ROOT / "03_coexpression_null/coexpression_summary.tsv", sep="\t", index=False)
    pd.DataFrame(selection_rows).to_csv(ROOT / "03_coexpression_null/coexpression_tolerance_ledger.tsv", sep="\t", index=False)

    # Three prespecified frozen controls; no new signature definitions.
    control_rows = []
    for sid in ["NEFTEL_AC_LIKE", "HALLMARK_EMT", "NEFTEL_NPC_LIKE"]:
        reference_genes = sig.loc[sig.signature_id.eq(sid), "gene"].astype(str).str.upper().drop_duplicates().tolist()
        genes = [g for g in reference_genes if g in z.columns]
        score = z[genes].mean(axis=1)
        pairs, caliper = match_ids(ng_ids, tr_ids, score)
        pm = pair_metrics(pairs, z, genes, score)
        obs_corr = float(pm.contribution_correlation.mean())
        lo, hi = bootstrap_mean(pm.contribution_correlation.to_numpy(float), rng)
        perm = label_permutation(strict_ids, len(ng_ids), score, z, genes, rng)
        control_rows.append({
            "signature_id": sid, "reference_gene_n": len(reference_genes), "measurable_gene_n": len(genes),
            "matched_pair_n": len(pairs), "caliper": caliper, "contribution_correlation": obs_corr,
            "bootstrap_CI_low": lo, "bootstrap_CI_high": hi, "permutation_null_mean": np.mean(perm),
            "empirical_P_lower_similarity": empirical_lower(perm, obs_corr), "role": "SUPPLEMENTARY_POSITIVE_CONTROL",
        })
    pd.DataFrame(control_rows).to_csv(ROOT / "04_cross_null_comparison/positive_control_signatures.tsv", sep="\t", index=False)

    comparison = pd.DataFrame([
        {"null": "EXACT91_EXPRESSION_MATCHED", "effective_n": len(eligible), "MES95_percentile": divergence_percentile, "empirical_P": exact_p, "decision": exact_decision},
        {"null": "COEXPRESSION_AWARE_SECONDARY", "effective_n": effective_n, "MES95_percentile": coex_percentile, "empirical_P": coex_p, "decision": coex_decision},
        {"null": "V2_1_95_MEASURED_GENE_HISTORICAL", "effective_n": 1527, "MES95_percentile": 99.148657, "empirical_P": .009162, "decision": "HISTORICAL_COMPARISON_ONLY"},
    ])
    comparison.to_csv(ROOT / "04_cross_null_comparison/cross_null_comparison.tsv", sep="\t", index=False)

    # Required table mirrors.
    for source, name in [
        (ROOT / "01_measurable_mes95/mes95_measurable_91.tsv", "mes95_measurable_91.tsv"),
        (ROOT / "01_measurable_mes95/mes95_missing_4.tsv", "mes95_missing_4.tsv"),
        (ROOT / "02_exact91_random_null/exact91_random_null.tsv", "exact91_random_null.tsv"),
        (ROOT / "02_exact91_random_null/exact91_random_summary.tsv", "exact91_random_summary.tsv"),
        (ROOT / "03_coexpression_null/coexpression_null.tsv", "coexpression_null.tsv"),
        (ROOT / "03_coexpression_null/coexpression_summary.tsv", "coexpression_summary.tsv"),
        (ROOT / "04_cross_null_comparison/positive_control_signatures.tsv", "positive_control_signatures.tsv"),
    ]:
        pd.read_csv(source, sep="\t").to_csv(ROOT / "06_tables" / name, sep="\t", index=False)

    receipt = {
        "seed": SEED, "bootstrap_n": NBOOT, "permutation_n": NPERM, "exact91_random_n": NRANDOM,
        "reference_gene_n": len(reference), "measurable_gene_n": len(measurable), "missing_genes": missing,
        "strict_RNA": {"NextGen": len(ng_ids), "traditional": len(tr_ids)},
        "MES95_matched_pairs": len(mes_pairs), "score_reconstruction_error": reconstruction_error,
        "exact91_decision": exact_decision, "coexpression_decision": coex_decision,
        "coexpression_tolerance": chosen, "coexpression_effective_n": effective_n,
        "dependency_branch": "NOT_TOUCHED_FROZEN",
    }
    (ROOT / "10_logs/phase1_null_validation_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
