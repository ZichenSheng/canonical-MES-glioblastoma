#!/usr/bin/env python3
"""Patient-level canonical MES compartment analysis from GBmap raw counts.

Cells are descriptive observations. Every inferential record is keyed by the
original study and donor. Neftel2019 is retained as discovery/descriptive and is
excluded from the deduplicated external-study summary.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path

DATA_ROOT = Path(os.environ.get("GBM_MES_DATA_ROOT", "data"))
RESULT_ROOT = Path(os.environ.get("GBM_MES_OUTPUT_ROOT", "results"))
REPO_ROOT = Path(os.environ.get("GBM_MES_REPO_ROOT", Path(__file__).resolve().parents[2]))

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import chi2, norm


GBMAP = DATA_ROOT / "single_cell" / "gbmap.h5ad"
OBS = DATA_ROOT / "prepared" / "gbmap" / "cell_metadata.tsv"
VAR = DATA_ROOT / "prepared" / "gbmap" / "gene_metadata.tsv"
SIGNATURES = REPO_ROOT / "resources" / "signatures" / "signature_gene_sets.tsv"
OUT = RESULT_ROOT / "single_cell" / "gbmap"
SEED = 2026073115
MIN_CELLS = 100
MIN_HIGH = 20
BLOCK = 5000
GBMAP_SHA256 = "459444da84efc45565ac9ccb68c4873dac033c07abfe395bc814b3a33f60c56c"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def compartment(cell_type: str, level4: str) -> str:
    ct = str(cell_type)
    l4 = str(level4).lower()
    if ct == "malignant cell":
        return "malignant"
    if ct in {"macrophage", "microglial cell", "monocyte", "dendritic cell"}:
        return "myeloid"
    if ct in {"mature T cell", "B cell", "plasma cell", "natural killer cell", "mast cell"}:
        return "lymphoid"
    if ct in {"endothelial cell", "mural cell"} or "vlmc" in l4 or "fibroblast" in l4:
        return "vascular_stromal"
    if ct in {"oligodendrocyte", "oligodendrocyte precursor cell"}:
        return "oligodendroglial"
    return "other_brain_or_unresolved"


def random_effect(y: np.ndarray, v: np.ndarray) -> dict[str, float]:
    ok = np.isfinite(y) & np.isfinite(v) & (v > 0)
    y, v = y[ok], v[ok]
    k = len(y)
    if k == 0:
        return {k0: np.nan for k0 in ["estimate", "se", "ci_low", "ci_high", "prediction_low", "prediction_high", "tau2", "I2", "Q", "Q_p"]} | {"k": 0}
    w = 1.0 / v
    fixed = float(np.sum(w * y) / np.sum(w))
    q = float(np.sum(w * (y - fixed) ** 2))
    c = float(np.sum(w) - np.sum(w**2) / np.sum(w))
    tau2 = max(0.0, (q - (k - 1)) / c) if k > 1 and c > 0 else 0.0
    wr = 1.0 / (v + tau2)
    mu = float(np.sum(wr * y) / np.sum(wr))
    se = float(math.sqrt(1.0 / np.sum(wr)))
    z = float(norm.ppf(0.975))
    pred = z * math.sqrt(tau2 + se**2) if k >= 3 else np.nan
    i2 = max(0.0, (q - (k - 1)) / q) * 100.0 if q > 0 and k > 1 else 0.0
    q_p = float(chi2.sf(q, k - 1)) if k > 1 else np.nan
    return {
        "k": k, "estimate": mu, "se": se, "ci_low": mu - z * se,
        "ci_high": mu + z * se, "prediction_low": mu - pred if np.isfinite(pred) else np.nan,
        "prediction_high": mu + pred if np.isfinite(pred) else np.nan,
        "tau2": tau2, "I2": i2, "Q": q, "Q_p": q_p,
    }


def decode_feature_names(var: pd.DataFrame) -> np.ndarray:
    for col in ("feature_name", "gene_symbol", "raw_var_id"):
        if col in var.columns:
            values = var[col].fillna("").astype(str).str.upper().str.strip().to_numpy()
            if np.count_nonzero(values != ""):
                return values
    raise RuntimeError("No usable raw feature-name field")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    script = Path(__file__).resolve()
    if sha256(GBMAP) != GBMAP_SHA256:
        raise SystemExit("GBmap SHA256 mismatch")

    obs = pd.read_csv(
        OBS, sep="\t", low_memory=False,
        usecols=["author", "donor_id", "disease", "cell_type", "annotation_level_4", "suspension_type", "assay"],
    )
    var = pd.read_csv(VAR, sep="\t", low_memory=False)
    if len(obs) != 338_564 or len(var) != 27_632:
        raise SystemExit("P0_STOP GBmap dimension mismatch")
    obs["valid"] = obs["disease"].eq("glioblastoma") & obs["donor_id"].notna() & obs["author"].notna()
    obs["study"] = obs["author"].astype(str)
    obs["patient_id"] = obs["study"] + "::" + obs["donor_id"].astype(str)
    obs["compartment"] = [compartment(a, b) for a, b in zip(obs["cell_type"], obs["annotation_level_4"])]
    obs["platform_class"] = np.where(obs["suspension_type"].astype(str).str.lower().eq("nucleus"), "snRNA", "scRNA")

    reg = pd.read_csv(SIGNATURES, sep="\t", low_memory=False)
    canonical = sorted(set(reg.loc[reg["signature_id"].eq("NEFTEL_MES_LIKE"), "gene_symbol_clean"].dropna().astype(str).str.upper().str.strip()))
    feature_names = decode_feature_names(var)
    gene_to_columns: dict[str, list[int]] = {}
    for index, gene in enumerate(feature_names):
        if gene in canonical:
            gene_to_columns.setdefault(gene, []).append(index)
    genes = sorted(gene_to_columns)
    if len(genes) < 20:
        raise SystemExit(f"P0_STOP canonical MES coverage too low: {len(genes)}")
    selected_columns = np.array([c for gene in genes for c in gene_to_columns[gene]], dtype=np.int64)
    col_to_gene = np.array([i for i, gene in enumerate(genes) for _ in gene_to_columns[gene]], dtype=np.int64)
    collapse = sparse.csr_matrix((np.ones(len(col_to_gene)), (np.arange(len(col_to_gene)), col_to_gene)), shape=(len(col_to_gene), len(genes)))

    n_cells = len(obs)
    normalized = np.full((n_cells, len(genes)), np.nan, dtype=np.float32)
    raw_relevant = np.zeros((n_cells, len(genes)), dtype=np.float32)
    library_size = np.zeros(n_cells, dtype=np.float64)
    with h5py.File(GBMAP, "r") as h5:
        raw = h5["raw/X"]
        indptr = raw["indptr"][:]
        shape = tuple(raw.attrs["shape"])
        if shape != (n_cells, len(var)):
            raise SystemExit(f"P0_STOP raw shape mismatch: {shape}")
        for start in range(0, n_cells, BLOCK):
            end = min(n_cells, start + BLOCK)
            data_start, data_end = int(indptr[start]), int(indptr[end])
            block = sparse.csr_matrix(
                (raw["data"][data_start:data_end], raw["indices"][data_start:data_end], indptr[start:end + 1] - data_start),
                shape=(end - start, shape[1]),
            )
            totals = np.asarray(block.sum(axis=1)).ravel().astype(np.float64)
            rel = (block[:, selected_columns] @ collapse).toarray().astype(np.float32)
            library_size[start:end] = totals
            raw_relevant[start:end] = rel
            scale = np.divide(10_000.0, totals, out=np.zeros_like(totals), where=totals > 0)
            normalized[start:end] = np.log1p(rel * scale[:, None]).astype(np.float32)
            print(f"GBMAP_BLOCK {start}:{end}", flush=True)

    # Gene-wise z scores are calculated within original study to avoid a direct
    # cross-platform scale comparison. High-score thresholds remain within donor.
    score = np.full(n_cells, np.nan, dtype=np.float64)
    for study, idx in obs.index[obs["valid"]].to_series().groupby(obs.loc[obs["valid"], "study"]).groups.items():
        ii = np.asarray(list(idx), dtype=np.int64)
        x = normalized[ii].astype(np.float64)
        mean = np.nanmean(x, axis=0)
        sd = np.nanstd(x, axis=0, ddof=1)
        keep = np.isfinite(sd) & (sd > 0)
        score[ii] = np.nanmean((x[:, keep] - mean[keep]) / sd[keep], axis=1)
    obs["canonical_MES"] = score

    valid = obs["valid"] & np.isfinite(obs["canonical_MES"])
    patients = obs.loc[valid, ["study", "patient_id", "platform_class"]].copy()
    patient_keys = sorted(patients["patient_id"].unique())
    patient_index = {p: i for i, p in enumerate(patient_keys)}
    gene_sums = np.zeros((len(patient_keys), len(genes)), dtype=np.float64)
    total_sums = np.zeros(len(patient_keys), dtype=np.float64)
    for p, frame in obs.loc[valid].groupby("patient_id", sort=True):
        ii = frame.index.to_numpy(dtype=np.int64)
        j = patient_index[p]
        gene_sums[j] = raw_relevant[ii].sum(axis=0)
        total_sums[j] = library_size[ii].sum()
    patient_log = np.log1p(gene_sums * np.divide(1e6, total_sums, out=np.zeros_like(total_sums), where=total_sums > 0)[:, None])
    patient_score = np.full(len(patient_keys), np.nan)
    patient_study = np.array([p.split("::", 1)[0] for p in patient_keys])
    for study in sorted(set(patient_study)):
        ii = np.where(patient_study == study)[0]
        x = patient_log[ii]
        sd = np.std(x, axis=0, ddof=1) if len(ii) > 1 else np.zeros(len(genes))
        keep = np.isfinite(sd) & (sd > 0)
        if np.any(keep):
            patient_score[ii] = np.mean((x[:, keep] - np.mean(x[:, keep], axis=0)) / sd[keep], axis=1)

    compartment_rows: list[dict] = []
    pseudo_rows: list[dict] = []
    all_compartments = ["malignant", "myeloid", "vascular_stromal", "lymphoid", "oligodendroglial", "other_brain_or_unresolved"]
    for p, frame in obs.loc[valid].groupby("patient_id", sort=True):
        n = len(frame)
        threshold = float(frame["canonical_MES"].quantile(0.75))
        high = frame["canonical_MES"].to_numpy() >= threshold
        n_high = int(high.sum())
        eligible = n >= MIN_CELLS and n_high >= MIN_HIGH
        platforms = sorted(frame["platform_class"].unique())
        platform = platforms[0] if len(platforms) == 1 else "mixed"
        study = str(frame["study"].iloc[0])
        pseudo_rows.append({
            "study": study, "patient_id": p, "platform_class": platform, "n_cells": n,
            "n_high": n_high, "high_threshold": threshold, "eligible": eligible,
            "canonical_MES_cell_median": float(frame["canonical_MES"].median()),
            "canonical_MES_pseudobulk": float(patient_score[patient_index[p]]), "gene_coverage": len(genes),
            "inference_unit": "patient", "discovery_role": "DISCOVERY_DESCRIPTIVE" if study == "Neftel2019" else "EXTERNAL_OR_REFERENCE_SOURCE",
        })
        for comp in all_compartments:
            is_comp = frame["compartment"].eq(comp).to_numpy()
            a = int(np.sum(high & is_comp)); b = int(np.sum(high & ~is_comp))
            c = int(np.sum(~high & is_comp)); d = int(np.sum(~high & ~is_comp))
            bg_fraction = float(np.mean(is_comp))
            hi_fraction = float(np.mean(is_comp[high])) if n_high else np.nan
            oe = hi_fraction / bg_fraction if bg_fraction > 0 else np.nan
            aa, bb, cc, dd = map(float, (a, b, c, d))
            correction = any(x == 0 for x in (aa, bb, cc, dd))
            if correction:
                aa += 0.5; bb += 0.5; cc += 0.5; dd += 0.5
            log_or = math.log((aa * dd) / (bb * cc))
            se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
            compartment_rows.append({
                "study": study, "patient_id": p, "platform_class": platform, "score": "canonical_MES",
                "compartment": comp, "n_cells": n, "n_high": n_high, "high_threshold": threshold,
                "background_count": int(is_comp.sum()), "high_count": a, "background_fraction": bg_fraction,
                "high_fraction": hi_fraction, "observed_expected": oe, "log2_observed_expected": math.log2(oe) if oe > 0 else np.nan,
                "odds_ratio": math.exp(log_or), "log_odds_ratio": log_or, "log_odds_se": se,
                "haldane_correction": correction, "eligible": eligible,
                "inference_unit": "patient", "pooled_cell_inference": False,
                "discovery_role": "DISCOVERY_DESCRIPTIVE" if study == "Neftel2019" else "EXTERNAL_OR_REFERENCE_SOURCE",
            })

    detail = pd.DataFrame(compartment_rows)
    pseudo = pd.DataFrame(pseudo_rows)
    eligible = detail.loc[detail["eligible"]].copy()
    study_rows: list[dict] = []
    for (study, comp), frame in eligible.groupby(["study", "compartment"], sort=True):
        fit = random_effect(frame["log_odds_ratio"].to_numpy(), frame["log_odds_se"].to_numpy() ** 2)
        study_rows.append({"level": "STUDY", "study": study, "omitted_study": "NONE", "platform_class": ";".join(sorted(frame["platform_class"].unique())), "compartment": comp, **fit})
    study_meta = pd.DataFrame(study_rows)

    overall_rows: list[dict] = []
    external_studies = sorted(set(study_meta["study"]) - {"Neftel2019"})
    for comp, frame in study_meta[study_meta["study"].isin(external_studies)].groupby("compartment", sort=True):
        fit = random_effect(frame["estimate"].to_numpy(), frame["se"].to_numpy() ** 2)
        overall_rows.append({"level": "EXTERNAL_STUDY_META", "study": "ALL_EXTERNAL", "omitted_study": "NONE", "platform_class": "ALL", "compartment": comp, **fit})
        for omit in external_studies:
            sub = frame.loc[frame["study"].ne(omit)]
            fit_o = random_effect(sub["estimate"].to_numpy(), sub["se"].to_numpy() ** 2)
            overall_rows.append({"level": "LOSO_EXTERNAL_STUDY_META", "study": "ALL_EXTERNAL", "omitted_study": omit, "platform_class": "ALL", "compartment": comp, **fit_o})
    overall = pd.DataFrame(overall_rows)
    meta = pd.concat([study_meta, overall], ignore_index=True)

    prediction = overall.loc[overall["level"].eq("EXTERNAL_STUDY_META")].copy()
    prediction["prediction_interval_status"] = np.where(prediction["k"] >= 3, "ESTIMABLE", "NOT_STABLE_K_LT_3")
    detail.to_csv(OUT / "PATIENT_LEVEL_CELL_COMPARTMENT.tsv", sep="\t", index=False)
    detail.to_csv(OUT / "PATIENT_LEVEL_SCORE_ENRICHMENT.tsv", sep="\t", index=False)
    pseudo.to_csv(OUT / "PATIENT_LEVEL_PSEUDOBULK.tsv", sep="\t", index=False)
    meta.to_csv(OUT / "PATIENT_LEVEL_META.tsv", sep="\t", index=False)
    prediction.to_csv(OUT / "PATIENT_LEVEL_PREDICTION_INTERVAL.tsv", sep="\t", index=False)

    qa = pd.DataFrame([
        {"check": "GBmap_sha256", "status": "PASS", "detail": GBMAP_SHA256},
        {"check": "GBmap_dimensions", "status": "PASS", "detail": f"cells={n_cells}; genes={len(var)}"},
        {"check": "canonical_gene_coverage", "status": "PASS", "detail": f"detected={len(genes)}; frozen_total={len(canonical)}"},
        {"check": "patient_identity", "status": "PASS", "detail": f"author::donor keys={len(patient_keys)}"},
        {"check": "eligible_patient_count", "status": "PASS" if pseudo["eligible"].sum() >= 5 else "FAIL", "detail": str(int(pseudo["eligible"].sum()))},
        {"check": "cell_independence", "status": "PASS", "detail": "cells used only for within-patient composition; formal rows are patient/study"},
        {"check": "GSE131928_role_isolation", "status": "PASS", "detail": "Neftel2019 retained as discovery/descriptive and excluded from external study meta"},
        {"check": "GSE174554_compartment_endpoint", "status": "NOT_EVALUABLE_IN_GBMAP", "detail": "handled independently and never represented as GBmap"},
        {"check": "contemporary_MES_pause", "status": "PASS", "detail": "no MES-Hyp/MES-Ast marker extraction or repeated benchmark performed"},
    ])
    qa.to_csv(OUT / "PATIENT_LEVEL_SINGLE_CELL_QA.tsv", sep="\t", index=False)
    freeze = {
        "schema_version": "1.1.0", "status": "COMPLETE", "seed": SEED,
        "script": str(script), "script_sha256": sha256(script),
        "input": str(GBMAP), "input_sha256": GBMAP_SHA256,
        "obs_sha256": sha256(OBS), "var_sha256": sha256(VAR), "signature_registry_sha256": sha256(SIGNATURES),
        "score": "canonical_MES", "normalization": "log1p(CP10K), gene-wise z within original study",
        "high_score": "within-patient upper quartile", "minimum_cells": MIN_CELLS, "minimum_high_cells": MIN_HIGH,
        "patient_key": "author::donor", "external_meta_excludes": ["Neftel2019"],
        "candidate_program_used": False, "contemporary_MES_repeat_used": False,
    }
    with (OUT / "PATIENT_LEVEL_SINGLE_CELL_FREEZE.json").open("w") as handle:
        json.dump(freeze, handle, indent=2, ensure_ascii=False)
    with (OUT / "PATIENT_LEVEL_SINGLE_CELL_GATE.md").open("w") as handle:
        handle.write("# Patient-level canonical MES analysis\n\nGBmap raw counts were scored with patients as the inferential unit. Neftel2019 is isolated from the deduplicated external-study meta-analysis. Cells are not treated as biological replicates. GSE174554 remains an independent analysis branch.\n")
    print(f"PATIENT_LEVEL_CANONICAL_COMPLETE patients={len(patient_keys)} eligible={int(pseudo['eligible'].sum())} studies={pseudo['study'].nunique()}")


if __name__ == "__main__":
    main()
