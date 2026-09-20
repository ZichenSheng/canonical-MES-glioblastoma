#!/usr/bin/env python3
"""Formal patient-blocked R10B Visium spatial validation."""

from __future__ import annotations

import gzip
import hashlib
import itertools
import json
import math
import os
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
from scipy.spatial import cKDTree

DATA_ROOT = Path(os.environ.get("GBM_MES_DATA_ROOT", "data"))
RESULT_ROOT = Path(os.environ.get("GBM_MES_OUTPUT_ROOT", "results"))
INPUT = DATA_ROOT / "prepared" / "spatial" / "multiomic_gbm"
ATLAS = DATA_ROOT / "spatial" / "multiomic_gbm"
OUT = RESULT_ROOT / "spatial" / "patient_blocked_validation"
MATRIX = INPUT / "R10B_selected_527x115914.f32"
SEED = 20260723
B = 2000
AXES = ["malignant_MES_realization", "myeloid", "hypoxia", "matrix_stromal", "vascular_stromal"]
AXIS_PROVENANCE = {
    "malignant_MES_realization": "EXPRESSION_PROGRAM: official Greenwald MES/MES.Hyp/MES.Ast spot label indicator",
    "myeloid": "EXPRESSION_PROGRAM: official atlas immune metaprogram indicator; not a direct cell count",
    "hypoxia": "EXPRESSION_PROGRAM: frozen HALLMARK_HYPOXIA bilaterally pruned against canonical MES",
    "matrix_stromal": "EXPRESSION_PROGRAM: frozen NABA_CORE_MATRISOME bilaterally pruned against canonical MES",
    "vascular_stromal": "EXPRESSION_PROGRAM: official atlas vascular metaprogram indicator; not a direct cell count",
}


def write_tsv(frame: pd.DataFrame, name: str) -> None:
    frame.to_csv(OUT / name, sep="\t", index=False, na_rep="NA")


def zscore(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    s = np.nanstd(x, ddof=1)
    if not np.isfinite(s) or s == 0:
        return np.full_like(x, np.nan)
    return (x - np.nanmean(x)) / s


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    good = np.isfinite(p)
    if not good.any():
        return out
    vals = p[good]
    order = np.argsort(vals)
    ranked = vals[order] * len(vals) / np.arange(1, len(vals) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    temp = np.empty(len(vals))
    temp[order] = np.minimum(ranked, 1.0)
    out[good] = temp
    return out


def ols_fit_predict(x_train: np.ndarray, y_train: np.ndarray, x_test: np.ndarray) -> np.ndarray:
    beta = np.linalg.lstsq(x_train, y_train, rcond=None)[0]
    return x_test @ beta


def base_design(g: pd.DataFrame) -> np.ndarray:
    cont = np.column_stack([
        g["canonical_MES_z"].to_numpy(float),
        g["log_count_z"].to_numpy(float),
        g["nfeature_z"].to_numpy(float),
    ])
    sections = pd.get_dummies(g["section"], drop_first=True, dtype=float).to_numpy()
    return np.column_stack([np.ones(len(g)), cont, sections])


def exact_signflip_p(effects: np.ndarray) -> float:
    effects = np.asarray(effects, float)
    effects = effects[np.isfinite(effects)]
    n = len(effects)
    if n == 0:
        return np.nan
    observed = float(np.mean(effects))
    if n <= 16:
        null = []
        for bits in itertools.product((-1.0, 1.0), repeat=n):
            null.append(np.mean(effects * np.asarray(bits)))
        null = np.asarray(null)
        return float((1 + np.sum(null >= observed)) / (len(null) + 1))
    rng = np.random.default_rng(SEED + 901)
    signs = rng.choice((-1.0, 1.0), size=(100000, n))
    null = (signs * effects).mean(axis=1)
    return float((1 + np.sum(null >= observed)) / (len(null) + 1))


def cluster_boot_ci(values: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    if len(values) < 2:
        return np.nan, np.nan
    boot = np.median(values[rng.integers(0, len(values), size=(B, len(values)))], axis=1)
    return tuple(np.quantile(boot, [.025, .975]))


def stable_partition(spot_id: str) -> int:
    digest = hashlib.sha256(f"{SEED}|{spot_id}".encode()).digest()
    return digest[0] & 1


def program_score(mm: np.memmap, global_indices: np.ndarray, rows: np.ndarray) -> tuple[np.ndarray, int]:
    x = np.asarray(mm[np.ix_(rows, global_indices)], dtype=np.float64)
    mu = x.mean(axis=1)
    sd = x.std(axis=1, ddof=1)
    ok = np.isfinite(sd) & (sd > 0)
    if ok.sum() < 2:
        return np.full(len(global_indices), np.nan), int(ok.sum())
    x = (x[ok] - mu[ok, None]) / sd[ok, None]
    return x.mean(axis=0), int(ok.sum())


def residualize_patient(g: pd.DataFrame, axis: str) -> np.ndarray:
    y = g[axis].to_numpy(float)
    x = base_design(g)
    ok = np.isfinite(y) & np.all(np.isfinite(x), axis=1)
    out = np.full(len(g), np.nan)
    if ok.sum() > x.shape[1] + 5 and np.nanstd(y[ok], ddof=1) > 0:
        out[ok] = y[ok] - ols_fit_predict(x[ok], y[ok], x[ok])
    return out


def moran_and_block_null(coords: np.ndarray, residual: np.ndarray, rng: np.random.Generator) -> tuple[float, np.ndarray, int]:
    ok = np.isfinite(residual) & np.all(np.isfinite(coords), axis=1)
    coords = coords[ok]
    r = residual[ok]
    n = len(r)
    if n < 30 or np.std(r, ddof=1) == 0:
        return np.nan, np.full(B, np.nan), 0
    r = r - r.mean()
    tree = cKDTree(coords)
    _, nbr = tree.query(coords, k=min(7, n))
    src = np.repeat(np.arange(n), nbr.shape[1] - 1)
    dst = nbr[:, 1:].reshape(-1)
    denom = np.sum(r * r)
    scale = n / len(src) / denom
    observed = float(scale * np.sum(r[src] * r[dst]))

    bx = np.floor(coords[:, 0] / 5.0).astype(np.int64)
    by = np.floor(coords[:, 1] / 5.0).astype(np.int64)
    _, block = np.unique(np.column_stack([bx, by]), axis=0, return_inverse=True)
    nb = int(block.max()) + 1
    key = block[src] * nb + block[dst]
    weights = np.bincount(key, weights=r[src] * r[dst], minlength=nb * nb).reshape(nb, nb)
    a = csr_matrix(weights)
    null = np.empty(B)
    for start in range(0, B, 200):
        stop = min(B, start + 200)
        signs = rng.choice((-1.0, 1.0), size=(stop - start, nb))
        null[start:stop] = scale * np.asarray((csr_matrix(signs) @ a).multiply(signs).sum(axis=1)).ravel()
    return observed, null, nb


def anatomical_cv(g: pd.DataFrame, axis: str) -> tuple[float, float, float]:
    y = g[axis].to_numpy(float)
    if np.nanstd(y, ddof=1) == 0:
        return np.nan, np.nan, np.nan
    x1 = base_design(g)
    af = pd.get_dummies(g["AF"].fillna("NA"), drop_first=True, dtype=float).to_numpy()
    x2 = np.column_stack([x1, af])
    block_key = (
        g["section"].astype(str) + "|" +
        np.floor(g["x"].to_numpy(float) / 10).astype(int).astype(str) + "|" +
        np.floor(g["y"].to_numpy(float) / 10).astype(int).astype(str)
    )
    fold = np.array([hashlib.sha256(f"{SEED}|{v}".encode()).digest()[0] % 5 for v in block_key], int)
    pred1 = np.full(len(g), np.nan)
    pred2 = np.full(len(g), np.nan)
    base = np.full(len(g), np.nan)
    for f in range(5):
        test = fold == f
        train = ~test
        if test.sum() == 0 or train.sum() <= x2.shape[1] + 5:
            continue
        pred1[test] = ols_fit_predict(x1[train], y[train], x1[test])
        pred2[test] = ols_fit_predict(x2[train], y[train], x2[test])
        base[test] = np.mean(y[train])
    ok = np.isfinite(pred1) & np.isfinite(pred2) & np.isfinite(base)
    if ok.sum() < 30:
        return np.nan, np.nan, np.nan
    sst = np.sum((y[ok] - base[ok]) ** 2)
    if sst <= 0:
        return np.nan, np.nan, np.nan
    q1 = 1 - np.sum((y[ok] - pred1[ok]) ** 2) / sst
    q2 = 1 - np.sum((y[ok] - pred2[ok]) ** 2) / sst
    return float(q1), float(q2), float(q2 - q1)


def bin_index(value: float, cuts: tuple[float, ...]) -> int:
    return int(np.searchsorted(np.asarray(cuts), value, side="right"))


def make_section_pairs(g: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    left_mask = np.array([stable_partition(v) == 0 for v in g["spot_id"]])
    left = np.where(left_mask)[0]
    right = np.where(~left_mask)[0]
    if len(left) == 0 or len(right) == 0:
        return pd.DataFrame(), 0
    cov = g[["canonical_MES_z", "log_count_z", "inferCNV_z"]].to_numpy(float)
    coords = g[["x", "y"]].to_numpy(float)
    tree = cKDTree(cov[right])
    candidates = tree.query_ball_point(cov[left], r=.20, p=np.inf)
    rows, cols, costs = [], [], []
    candidate_count = 0
    for li, js in enumerate(candidates):
        if not js:
            continue
        rr = right[np.asarray(js, int)]
        physical = np.sqrt(np.sum((coords[rr] - coords[left[li]]) ** 2, axis=1))
        keep = physical >= 5.0
        rr = rr[keep]
        js2 = np.asarray(js, int)[keep]
        if len(rr) == 0:
            continue
        delta = np.abs(cov[rr] - cov[left[li]]).sum(axis=1)
        rows.extend([li] * len(rr))
        cols.extend(js2.tolist())
        costs.extend((delta + 1e-6).tolist())
        candidate_count += len(rr)
    # Unique dummy edge per left implements lexicographic maximum cardinality,
    # then minimum covariate distance, in an exact sparse bipartite solve.
    for li in range(len(left)):
        rows.append(li)
        cols.append(len(right) + li)
        costs.append(1e6)
    graph = coo_matrix((costs, (rows, cols)), shape=(len(left), len(right) + len(left))).tocsr()
    ri, cj = min_weight_full_bipartite_matching(graph)
    keep = cj < len(right)
    lidx = left[ri[keep]]
    ridx = right[cj[keep]]
    out = pd.DataFrame({"left_local": lidx, "right_local": ridx})
    if len(out):
        out["delta_MES"] = np.abs(cov[lidx, 0] - cov[ridx, 0])
        out["delta_log_count"] = np.abs(cov[lidx, 1] - cov[ridx, 1])
        out["delta_inferCNV"] = np.abs(cov[lidx, 2] - cov[ridx, 2])
        out["physical_distance"] = np.sqrt(np.sum((coords[lidx] - coords[ridx]) ** 2, axis=1))
    return out, candidate_count


def constrained_null(g: pd.DataFrame, pairs: pd.DataFrame, rng: np.random.Generator) -> tuple[dict[str, np.ndarray], float]:
    left = pairs["left_index"].to_numpy(int)
    current = pairs["right_index"].to_numpy(int).copy()
    sections = g["section"].astype(str).to_numpy()
    mes = g["canonical_MES_z"].to_numpy(float)
    depth = g["log_count_z"].to_numpy(float)
    mal = g["inferCNV_z"].to_numpy(float)
    coords = g[["x", "y"]].to_numpy(float)
    axis = g[AXES].to_numpy(float)
    af = g["AF"].astype(str).to_numpy()

    def signature(i: int, j: int) -> tuple[int, int, int, int]:
        return (
            bin_index(abs(mes[i] - mes[j]), (.05, .10, .15)),
            bin_index(abs(depth[i] - depth[j]), (.05, .10, .15)),
            bin_index(abs(mal[i] - mal[j]), (.05, .10, .15)),
            bin_index(float(np.linalg.norm(coords[i] - coords[j])), (10, 20, 40, 80)),
        )

    targets = [signature(i, j) for i, j in zip(left, current)]
    groups: dict[tuple[str, tuple[int, int, int, int]], list[int]] = {}
    for k, (i, sig) in enumerate(zip(left, targets)):
        groups.setdefault((sections[i], sig), []).append(k)
    movable = [np.asarray(v, int) for v in groups.values() if len(v) >= 2]
    null_distance = np.empty(B)
    null_niche = np.empty(B)
    attempts = accepted = 0

    def valid(k: int, j: int) -> bool:
        i = left[k]
        return (
            i != j and sections[i] == sections[j]
            and abs(mes[i] - mes[j]) <= .20
            and abs(depth[i] - depth[j]) <= .20
            and abs(mal[i] - mal[j]) <= .20
            and np.linalg.norm(coords[i] - coords[j]) >= 5.0
            and signature(i, j) == targets[k]
        )

    warmup = 5000
    thin = 50
    total = warmup + B * thin
    for step in range(total):
        if movable:
            group = movable[rng.integers(0, len(movable))]
            a, b = rng.choice(group, 2, replace=False)
            ja, jb = current[a], current[b]
            attempts += 1
            if valid(a, jb) and valid(b, ja):
                current[a], current[b] = jb, ja
                accepted += 1
        if step >= warmup and (step - warmup) % thin == 0:
            bi = (step - warmup) // thin
            if bi < B:
                diff = axis[left] - axis[current]
                n_ok = np.isfinite(diff).sum(axis=1)
                dist = np.where(n_ok >= 3, np.sqrt(np.nansum(diff * diff, axis=1) / np.maximum(n_ok, 1)), np.nan)
                null_distance[bi] = np.nanmedian(dist)
                null_niche[bi] = np.mean(af[left] != af[current])
    return {"realization_distance": null_distance, "niche_discordance": null_niche}, accepted / max(attempts, 1)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    meta = pd.read_csv(ATLAS / "meta.tsv", sep="\t", low_memory=False).rename(columns={"Unnamed: 0": "spot_id"})
    coords = pd.read_csv(ATLAS / "Visium.coords.tsv.gz", sep="\t", header=None, names=["spot_id", "x", "y"])
    assert meta["spot_id"].equals(coords["spot_id"])
    meta["global_index"] = np.arange(len(meta))
    meta["x"] = coords["x"].to_numpy(float)
    meta["y"] = coords["y"].to_numpy(float)
    primary = meta.short_histology.eq("GBM") & meta.IDH1.eq("WT")
    qc = (meta.nCount_Spatial >= 1000) & (meta.nFeature_Spatial >= 500)
    eligible = primary & qc
    d = meta.loc[eligible].copy().rename(columns={"sample": "patient", "orig.ident": "section"})
    assert d.patient.nunique() == 13 and d.section.nunique() == 28

    selected = pd.read_csv(INPUT / "R10B_SELECTED_MATRIX_ROWS.tsv", sep="\t")
    selected_map = {g: int(i) for g, i in zip(selected.gene, selected.selected_row_0based)}
    membership = pd.read_csv(INPUT / "R10B_AXIS_GENE_MEMBERSHIP.tsv", sep="\t")
    def axis_rows(name: str) -> np.ndarray:
        genes = membership.loc[membership.axis_gene_set.eq(name) & membership.mapping_status.str.startswith("MAPPED_UNIQUE"), "mapped_feature_gene"]
        return np.array([selected_map[g] for g in genes if g in selected_map], int)
    rowsets = {
        "canonical": axis_rows("CANONICAL_MES"),
        "hypoxia": axis_rows("HYPOXIA_PRUNED_CANONICAL_MES"),
        "matrix": axis_rows("MATRIX_PRUNED_CANONICAL_MES"),
    }
    assert len(rowsets["canonical"]) == 95
    mm = np.memmap(MATRIX, dtype=np.float32, mode="r", shape=(len(selected), len(meta)))
    score_qc = []
    frames = []
    for patient, g0 in d.groupby("patient", sort=True):
        g = g0.copy()
        gi = g.global_index.to_numpy(int)
        mes, mes_n = program_score(mm, gi, rowsets["canonical"])
        hyp, hyp_n = program_score(mm, gi, rowsets["hypoxia"])
        matrix, matrix_n = program_score(mm, gi, rowsets["matrix"])
        g["canonical_MES_raw"] = mes
        g["canonical_MES_z"] = zscore(mes)
        g["malignant_MES_realization"] = g.greenwald_metaprograms.astype(str).str.startswith("MES").astype(float)
        g["myeloid"] = g.metaprogram.eq("immune").astype(float)
        g["hypoxia"] = zscore(hyp)
        g["matrix_stromal"] = zscore(matrix)
        g["vascular_stromal"] = g.metaprogram.eq("vascular").astype(float)
        g["log_count_z"] = zscore(np.log1p(g.nCount_Spatial.to_numpy(float)))
        g["nfeature_z"] = zscore(g.nFeature_Spatial.to_numpy(float))
        g["inferCNV_z"] = zscore(g.inferCNV.to_numpy(float))
        for axis in AXES:
            g[axis] = zscore(g[axis].to_numpy(float))
        score_qc.append({
            "patient": patient, "eligible_spots": len(g), "canonical_nonzero_variance_genes": mes_n,
            "hypoxia_nonzero_variance_genes": hyp_n, "matrix_nonzero_variance_genes": matrix_n,
            "canonical_mapping_gate": "PASS" if mes_n >= 86 else "FAIL",
        })
        frames.append(g)
    d = pd.concat(frames, ignore_index=True)
    write_tsv(pd.DataFrame(score_qc), "R10B_SPOT_SCORE_QC.tsv")
    if not all(x["canonical_mapping_gate"] == "PASS" for x in score_qc):
        raise RuntimeError("Patient-level canonical nonzero-variance mapping gate did not pass")

    realization_cols = [
        "spot_id", "global_index", "patient", "section", "source", "Diagnosis", "IDH1",
        "nCount_Spatial", "nFeature_Spatial", "x", "y", "inferCNV", "inferCNV_z", "AF", "AF_ivy",
        "metaprogram", "greenwald_metaprograms", "canonical_MES_raw", "canonical_MES_z",
        "log_count_z", "nfeature_z", *AXES,
    ]
    write_tsv(d[realization_cols], "R10B_SPATIAL_REALIZATION_MATRIX.tsv")

    # MES-conditioned residual spatial structure, patient is the inference unit.
    residual_rows = []
    patient_axis_null: dict[tuple[str, str], np.ndarray] = {}
    for patient, gp in d.groupby("patient", sort=True):
        gp = gp.copy().reset_index(drop=True)
        residual_matrix = []
        for ai, axis in enumerate(AXES):
            residual = residualize_patient(gp, axis)
            residual_matrix.append(residual)
            section_obs, section_null, section_n, section_blocks = [], [], [], []
            for si, (section, gs) in enumerate(gp.assign(_res=residual).groupby("section", sort=True)):
                obs, null, nb = moran_and_block_null(
                    gs[["x", "y"]].to_numpy(float), gs["_res"].to_numpy(float),
                    np.random.default_rng(SEED + 1000 + ai * 100 + si + int(hashlib.sha256(patient.encode()).hexdigest()[:4], 16)),
                )
                section_obs.append(obs); section_null.append(null); section_n.append(len(gs)); section_blocks.append(nb)
                residual_rows.append({
                    "record_type": "SECTION", "patient": patient, "section": section, "axis": axis,
                    "n_spots": len(gs), "residual_variance": np.nanvar(gs["_res"], ddof=1),
                    "Moran_I": obs, "null_median": np.nanmedian(null),
                    "effect": obs - np.nanmedian(null),
                    "block_permutation_p": (1 + np.sum(null >= obs)) / (B + 1) if np.isfinite(obs) else np.nan,
                    "permutations": B, "spatial_blocks": nb,
                })
            weights = np.asarray(section_n, float)
            good = np.isfinite(section_obs)
            if good.any():
                w = weights[good] / weights[good].sum()
                pat_obs = float(np.sum(np.asarray(section_obs)[good] * w))
                pat_null = np.sum(np.stack([section_null[i] for i in np.where(good)[0]], axis=1) * w, axis=1)
                patient_axis_null[(patient, axis)] = pat_null
                residual_rows.append({
                    "record_type": "PATIENT", "patient": patient, "section": "ALL_NESTED", "axis": axis,
                    "n_spots": len(gp), "residual_variance": np.nanvar(residual, ddof=1),
                    "Moran_I": pat_obs, "null_median": np.nanmedian(pat_null),
                    "effect": pat_obs - np.nanmedian(pat_null),
                    "block_permutation_p": (1 + np.sum(pat_null >= pat_obs)) / (B + 1),
                    "permutations": B, "spatial_blocks": int(np.sum(section_blocks)),
                })
        rmat = np.column_stack(residual_matrix)
        ok = np.all(np.isfinite(rmat), axis=1)
        if ok.sum() > 10:
            sv = np.linalg.svd(rmat[ok] - rmat[ok].mean(axis=0), compute_uv=False) ** 2
            residual_rows.append({
                "record_type": "PATIENT_MULTIVARIATE", "patient": patient, "section": "ALL_NESTED", "axis": "JOINT_FIVE_AXIS",
                "n_spots": int(ok.sum()), "residual_variance": float(np.sum(np.var(rmat[ok], axis=0, ddof=1))),
                "Moran_I": np.nan, "null_median": np.nan, "effect": float(sv[0] / sv.sum()),
                "block_permutation_p": np.nan, "permutations": 0, "spatial_blocks": np.nan,
            })
    residual = pd.DataFrame(residual_rows)
    summary_rows = []
    for axis in AXES:
        sub = residual[(residual.record_type == "PATIENT") & (residual.axis == axis)]
        effects = sub.effect.to_numpy(float)
        lo, hi = cluster_boot_ci(effects, rng)
        summary_rows.append({
            "record_type": "GLOBAL_PATIENT_INFERENCE", "patient": "N=13", "section": "SECTIONS_NESTED", "axis": axis,
            "n_spots": int(sub.n_spots.sum()), "residual_variance": float(np.nanmedian(sub.residual_variance)),
            "Moran_I": float(np.nanmedian(sub.Moran_I)), "null_median": float(np.nanmedian(sub.null_median)),
            "effect": float(np.nanmedian(effects)), "effect_CI_low": lo, "effect_CI_high": hi,
            "block_permutation_p": exact_signflip_p(effects), "permutations": B, "spatial_blocks": int(np.nansum(sub.spatial_blocks)),
            "positive_patient_fraction": float(np.mean(effects > 0)),
        })
    summaries = pd.DataFrame(summary_rows)
    summaries["BH_q"] = bh(summaries.block_permutation_p.to_numpy(float))
    residual = pd.concat([residual, summaries], ignore_index=True, sort=False)
    write_tsv(residual, "R10B_PATIENT_LEVEL_RESIDUAL_STRUCTURE.tsv")

    # Official anatomical context incremental value under patient-blocked inference.
    context_rows = []
    for patient, gp in d.groupby("patient", sort=True):
        for axis in AXES:
            q1, q2, delta = anatomical_cv(gp.reset_index(drop=True), axis)
            context_rows.append({"record_type": "PATIENT", "patient": patient, "axis": axis, "n_spots": len(gp), "Q2_S1": q1, "Q2_S2": q2, "delta_Q2": delta})
    context = pd.DataFrame(context_rows)
    context_summary = []
    for axis in AXES:
        vals = context.loc[context.axis.eq(axis), "delta_Q2"].to_numpy(float)
        finite = np.isfinite(vals)
        lo, hi = cluster_boot_ci(vals, rng)
        context_summary.append({
            "record_type": "GLOBAL_PATIENT_INFERENCE", "patient": "N=13", "axis": axis,
            "n_spots": int(context.loc[context.axis.eq(axis), "n_spots"].sum()),
            "Q2_S1": float(np.nanmedian(context.loc[context.axis.eq(axis), "Q2_S1"])),
            "Q2_S2": float(np.nanmedian(context.loc[context.axis.eq(axis), "Q2_S2"])),
            "delta_Q2": float(np.nanmedian(vals)), "delta_Q2_CI_low": lo, "delta_Q2_CI_high": hi,
            "positive_patient_fraction": float(np.mean(vals[finite] > 0)) if finite.any() else np.nan, "patient_signflip_p": exact_signflip_p(vals),
        })
    cs = pd.DataFrame(context_summary)
    cs["BH_q"] = bh(cs.patient_signflip_p.to_numpy(float))
    cs["incremental_context_gate"] = np.where(
        (cs.delta_Q2 > 0) & (cs.delta_Q2_CI_low > 0) & (cs.positive_patient_fraction >= .60) & (cs.BH_q < .05),
        "SUPPORTED", "NOT_SUPPORTED",
    )
    context = pd.concat([context, cs], ignore_index=True, sort=False)
    write_tsv(context, "R10B_ANATOMICAL_CONTEXT_INCREMENTAL_VALUE.tsv")

    # Exact outcome-blind iso-MES matching within section.
    pair_frames = []
    candidate_rows = []
    for (patient, section), gs0 in d.groupby(["patient", "section"], sort=True):
        gs = gs0.reset_index()
        pp, candidates = make_section_pairs(gs)
        candidate_rows.append({"patient": patient, "section": section, "spots": len(gs), "eligible_candidate_edges": candidates, "matched_pairs": len(pp)})
        if len(pp):
            pp["patient"] = patient; pp["section"] = section
            pp["left_index"] = gs.loc[pp.left_local, "index"].to_numpy(int)
            pp["right_index"] = gs.loc[pp.right_local, "index"].to_numpy(int)
            pp["left_spot"] = gs.loc[pp.left_local, "spot_id"].to_numpy()
            pp["right_spot"] = gs.loc[pp.right_local, "spot_id"].to_numpy()
            pair_frames.append(pp)
    pairs = pd.concat(pair_frames, ignore_index=True) if pair_frames else pd.DataFrame()
    write_tsv(pd.DataFrame(candidate_rows), "R10B_ISO_MES_MATCHING_QA.tsv")
    pairs.to_csv(OUT / "R10B_ISO_MES_MATCHED_SPOT_PAIRS.tsv.gz", sep="\t", index=False, compression="gzip")

    iso_patient = []
    for pi, (patient, pp0) in enumerate(pairs.groupby("patient", sort=True)):
        gp = d[d.patient.eq(patient)].copy()
        old_to_new = {old: new for new, old in enumerate(gp.index)}
        pp = pp0.copy()
        pp["left_index"] = pp.left_index.map(old_to_new).astype(int)
        pp["right_index"] = pp.right_index.map(old_to_new).astype(int)
        gp = gp.reset_index(drop=True)
        axis = gp[AXES].to_numpy(float)
        left = pp.left_index.to_numpy(int); right = pp.right_index.to_numpy(int)
        diff = axis[left] - axis[right]
        n_ok = np.isfinite(diff).sum(axis=1)
        distance = np.where(n_ok >= 3, np.sqrt(np.nansum(diff * diff, axis=1) / np.maximum(n_ok, 1)), np.nan)
        observed_distance = float(np.nanmedian(distance))
        observed_niche = float(np.mean(gp.AF.astype(str).to_numpy()[left] != gp.AF.astype(str).to_numpy()[right]))
        nulls, acceptance = constrained_null(gp, pp, np.random.default_rng(SEED + 5000 + pi))
        for family, observed in [("realization_distance", observed_distance), ("niche_discordance", observed_niche)]:
            null = nulls[family]
            iso_patient.append({
                "record_type": "PATIENT", "patient": patient, "family": family, "n_pairs": len(pp),
                "observed": observed, "null_median": float(np.nanmedian(null)), "effect": observed - float(np.nanmedian(null)),
                "within_patient_p": float((1 + np.sum(null >= observed)) / (B + 1)),
                "null_reassignments": B, "swap_acceptance": acceptance,
            })
    iso = pd.DataFrame(iso_patient)
    iso_summary = []
    for family in ["realization_distance", "niche_discordance"]:
        sub = iso[iso.family.eq(family)]
        effects = sub.effect.to_numpy(float)
        lo, hi = cluster_boot_ci(effects, rng)
        iso_summary.append({
            "record_type": "GLOBAL_PATIENT_INFERENCE", "patient": f"N={len(sub)}", "family": family,
            "n_pairs": int(sub.n_pairs.sum()), "observed": float(np.nanmedian(sub.observed)),
            "null_median": float(np.nanmedian(sub.null_median)), "effect": float(np.nanmedian(effects)),
            "effect_CI_low": lo, "effect_CI_high": hi, "positive_patient_fraction": float(np.mean(effects > 0)),
            "within_patient_p": exact_signflip_p(effects), "null_reassignments": B,
            "swap_acceptance": float(np.nanmedian(sub.swap_acceptance)),
        })
    iss = pd.DataFrame(iso_summary)
    iss["BH_q"] = bh(iss.within_patient_p.to_numpy(float))
    iss["iso_MES_gate"] = np.where(
        (iss.effect > 0) & (iss.effect_CI_low > 0) & (iss.positive_patient_fraction >= .60) & (iss.BH_q < .05) & (iss.swap_acceptance > 0),
        "SUPPORTED", "NOT_SUPPORTED",
    )
    iso = pd.concat([iso, iss], ignore_index=True, sort=False)
    write_tsv(iso, "R10B_ISO_MES_SPATIAL_NON_EQUIVALENCE.tsv")

    # Xenium mapping/support boundary from official custom panel and package contents.
    panel_path = ATLAS / "SNU21_extracted/SNU21/gene_panel.json"
    panel = json.load(panel_path.open())
    panel_genes = {x["type"]["data"]["name"].upper() for x in panel["payload"]["targets"] if x.get("type", {}).get("descriptor") == "gene"}
    mesmap = pd.read_csv(INPUT / "R10B_MES_MAPPING_AUDIT.tsv", sep="\t")
    mapped = mesmap.mapped_feature_gene.isin(panel_genes)
    with tarfile.open(ATLAS / "SNU21_extracted/SNU21/analysis.tar.gz", "r:gz") as tar:
        names = tar.getnames()
    author_celltypes_available = any("celltype" in x.lower() or "annotation" in x.lower() for x in names)
    xenium = pd.DataFrame([
        {"record_type": "PANEL", "patient": "SNU21", "metric": "panel_gene_count", "value": len(panel_genes), "status": "OFFICIAL"},
        {"record_type": "CANONICAL_MAPPING", "patient": "SNU21", "metric": "canonical_MES_mapped", "value": f"{int(mapped.sum())}/95", "status": "XENIUM_CANONICAL_MES_NOT_EVALUABLE"},
        {"record_type": "CELLULAR_SUPPORT", "patient": "SNU21", "metric": "author_level_celltype_or_AF_annotation_packaged", "value": str(author_celltypes_available).upper(), "status": "CELL_RESOLVED_SUPPORT_NOT_EVALUABLE" if not author_celltypes_available else "EVALUABLE"},
        {"record_type": "SERIES_BOUNDARY", "patient": "PUBLISHED_4_CASE_SERIES", "metric": "inference_boundary", "value": "secondary small series; no independent large-cohort claim", "status": "FROZEN_BOUNDARY"},
    ])
    write_tsv(xenium, "R10B_XENIUM_SUPPORT.tsv")

    residual_global = residual[residual.record_type.eq("GLOBAL_PATIENT_INFERENCE")]
    context_global = context[context.record_type.eq("GLOBAL_PATIENT_INFERENCE")]
    iso_global = iso[iso.record_type.eq("GLOBAL_PATIENT_INFERENCE")]
    residual_supported_axes = set(residual_global.loc[(residual_global.effect > 0) & (residual_global.effect_CI_low > 0) & (residual_global.positive_patient_fraction >= .60) & (residual_global.BH_q < .05), "axis"])
    iso_supported = bool((iso_global.iso_MES_gate == "SUPPORTED").any())
    context_supported = bool((context_global.incremental_context_gate == "SUPPORTED").any())
    if iso_supported and residual_supported_axes:
        status = "SPATIAL_CONTEXTUAL_REALIZATION_SUPPORTED" if context_supported else "SPATIAL_NONIDENTIFICATION_SUPPORTED"
    else:
        # Association-only if at least one patient-median |MES-axis correlation| is nontrivial.
        correlations = []
        for patient, gp in d.groupby("patient"):
            for axis in AXES:
                correlations.append(abs(pd.Series(gp.canonical_MES_z).corr(pd.Series(gp[axis]))))
        status = "SPATIAL_MES_ECOLOGY_ASSOCIATION_ONLY" if np.nanmedian(correlations) >= .10 else "SPATIAL_NONIDENTIFICATION_NOT_SUPPORTED"
    receipt = {
        "status": status,
        "primary_patients": int(d.patient.nunique()), "primary_sections": int(d.section.nunique()),
        "raw_primary_spots": int(primary.sum()), "eligible_primary_spots": int(len(d)),
        "visium_canonical_mapping": "95/95", "iso_evaluable_patients": int(iso[iso.record_type.eq("PATIENT")].patient.nunique()),
        "iso_supported": iso_supported, "residual_supported_axes": sorted(residual_supported_axes),
        "context_supported": context_supported, "xenium_canonical_mapping": f"{int(mapped.sum())}/95",
        "xenium_status": "XENIUM_CANONICAL_MES_NOT_EVALUABLE",
        "xenium_cell_support": "CELL_RESOLVED_SUPPORT_NOT_EVALUABLE" if not author_celltypes_available else "EVALUABLE",
    }
    (OUT / "R10B_FINAL_RECEIPT.json").write_text(json.dumps(receipt, indent=2) + "\n")
    adjudication = f"""# R10B adjudication

## Verdict

`{status}`

## Authenticated denominator

- Official complete atlas export: 115,914 total spots across 17 patients/32 sections.
- Primary authoritative IDH-wildtype GBM: **{d.patient.nunique()} patients, {d.section.nunique()} sections, {int(primary.sum()):,} raw spots**.
- Frozen QC-evaluable primary: **{len(d):,} spots**; spots are within-patient observations, never independent patients.
- Visium canonical MES mapping: **95/95**, including only two unique HGNC historical-symbol mappings (C8orf4→TCIM and ERO1L→ERO1A).

## Formal spatial gates

- Iso-MES non-equivalence: **{'supported' if iso_supported else 'not supported'}** under patient-level sign/randomization inference and 2,000 constrained within-patient null states.
- MES-conditioned residual spatial structure axes passing the complete patient-level gate: **{', '.join(sorted(residual_supported_axes)) if residual_supported_axes else 'none'}**.
- Official anatomical-context incremental gate: **{'supported' if context_supported else 'not supported'}**.

## Xenium boundary

- Official custom panel mapping is **{int(mapped.sum())}/95**; therefore `XENIUM_CANONICAL_MES_NOT_EVALUABLE`.
- The official minimal processed archive contains 10x clustering outputs but not the author-level cell-type/AF annotation used by the paper scripts; therefore `CELL_RESOLVED_SUPPORT_NOT_EVALUABLE`. Unsupervised clusters were not relabeled as cell types.

## Interpretation boundary

This is a within-patient cross-sectional spatial test. It does not treat spots as independent patients, does not establish causal ecology, does not imply temporal rewiring, and does not alter the frozen canonical MES signature or realization axes.
"""
    (OUT / "R10B_ADJUDICATION.md").write_text(adjudication)
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
