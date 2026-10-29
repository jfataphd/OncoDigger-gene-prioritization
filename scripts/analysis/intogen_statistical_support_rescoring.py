"""
intogen_statistical_support_rescoring.py -- sensitivity analysis: how much of the
IntOGen concordance survives when the reference is restricted to catalogue rows that
meet IntOGen's own non-knowledge-based statistical admission rule?

Why: IntOGen is not literature-orthogonal. Its catalogue admission uses the Cancer Gene
Census (CGC, literature-curated) to set method weights, to waive the two-method
requirement, to rescue genes at a relaxed FDR (Tier 2) and to exempt genes from two
filters, and a CancerMine snapshot as a removal-only gate. See Methods 2.5 and
Results 5.9 of the manuscript.

Reference definitions (all from the bundled release 2024.09.20 Compendium_Cancer_Genes.tsv):
  published   every catalogue row                                   (reproduces 177/190, 144/190)
  statistical catalogue rows with QVALUE_COMBINATION <= 0.05 AND >= 2 individual
              significant methods (the METHODS entry "combination" is the weighted
              Stouffer p-value combination, not an individual method, and is not counted).
              This is the rule IntOGen applies to a gene that is NOT in the CGC, so a
              row that passes it needs no CGC-specific admission route.

This is a sensitivity analysis, NOT an independent reference: IntOGen's method weights
are CGC-calibrated and cannot be removed without re-running IntOGen on the raw cohorts.

Per-cancer bootstrap: 2000 resamples of the 19 cancers, seed 42 (same scheme as
intogen_precision_bootstrap_ci.py).

Inputs  (all in this repository):
  data/canonical/precision_at_k_detail.csv, data/derived/*.csv,
  data/derived/intogen_rescoring/*_detail.csv,
  data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv
Outputs (data/derived/intogen_rescoring/):
  statistical_support_reference_summary.csv
  statistical_support_reference_paired_vs_oncodigger.csv
  statistical_support_reference_by_cancer_oncodigger.csv

Run: python -X utf8 scripts/analysis/intogen_statistical_support_rescoring.py
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DERIVED = ROOT / "data" / "derived"
RESC = DERIVED / "intogen_rescoring"
INTOGEN_TSV = ROOT / "data" / "external" / "intogen" / "drivers" / "2024-06-18_IntOGen-Drivers" / "Compendium_Cancer_Genes.tsv"
Q_THRESHOLD = 0.05
MIN_METHODS = 2
BOOT_B, BOOT_SEED = 2000, 42

# Supplementary Table S8 (published) -- identical to verify_headline_results.py.
CANCER_MAP = {
    "Breast Cancer": {"BRCA"},
    "Lung Cancer": {"LUAD", "LUSC", "NSCLC", "SCLC"},
    "Colorectal Cancer": {"COAD", "READ", "COADREAD"},
    "Prostate Cancer": {"PRAD"},
    "Melanoma": {"MEL", "SKCM", "UM"},
    "Bladder Cancer": {"BLCA", "UTUC"},
    "Kidney Cancer": {"CCRCC", "CHRCC", "PRCC", "RCC", "WT"},
    "Pancreatic Cancer": {"PAAD", "PANET"},
    "Liver Cancer": {"HCC", "CHOL", "LIHB"},
    "Stomach Cancer": {"STAD"},
    "Esophageal Cancer": {"ESCA", "ESCC"},
    "Ovarian Cancer": {"OVT"},
    "Endometrial Cancer": {"UCEC", "UCS"},
    "Cervical Cancer": {"CESC", "CEAD"},
    "Thyroid Cancer": {"WDTC"},
    "Brain Cancer": {"GB", "GBM", "HGGNOS", "LGGNOS", "PAST", "MBL", "EPM", "ATRT"},
    "Leukemia": {"ALL", "AML", "CLLSLL", "MDS"},
    "Lymphoma": {"BL", "DLBCLNOS", "NHL", "MLYM"},
    "Myeloma": {"PCM"},
}
CANCERS = list(CANCER_MAP)

# name -> (kind, file)
SYSTEMS = {
    "OncoDigger": ("canon", None),
    "Pooled cancer-blind control": ("pooled", None),
    "CancerMine": ("file", "cancermine_top10.csv"),
    "ChatGPT (one response)": ("file", "chatgpt_top10.csv"),
    "DISEASES 2.0": ("file", "diseases_top10.csv"),
    "PubTator 3.0": ("file", "pubtator_v2_top10.csv"),
    "Matched frequency-only (clean vocabulary)": ("resc", "matched_frequency_baseline_detail.csv"),
    "Frequency baseline (aliases)": ("file", "freq_baseline_top50_aliases.csv"),
    "Frequency baseline (primary symbols)": ("file", "freq_baseline_top50_primary.csv"),
}


def n_individual_methods(methods: str) -> int:
    return len([m for m in methods.split(",") if m and m != "combination"])


def load_references():
    rows = list(csv.DictReader(open(INTOGEN_TSV, encoding="utf-8-sig"), delimiter="\t"))
    refs = {}
    for name, keep in (
        ("published", lambda r: True),
        ("statistical", lambda r: float(r["QVALUE_COMBINATION"]) <= Q_THRESHOLD
         and n_individual_methods(r["METHODS"]) >= MIN_METHODS),
    ):
        sel = [r for r in rows if keep(r)]
        any_genes = {r["SYMBOL"].strip().upper() for r in sel}
        exact = set()
        for r in sel:
            for cancer, codes in CANCER_MAP.items():
                if r["CANCER_TYPE"].strip() in codes:
                    exact.add((r["SYMBOL"].strip().upper(), cancer))
        refs[name] = (any_genes, exact, len(sel))
    return refs


def load_system(kind: str, fname: str | None) -> pd.DataFrame:
    if kind == "canon":
        t = pd.read_csv(ROOT / "data" / "canonical" / "precision_at_k_detail.csv")
    elif kind == "pooled":
        g = pd.read_csv(DERIVED / "pooled_control_mutation.csv").head(10)
        t = pd.DataFrame([(c, r, x) for c in CANCERS for r, x in zip(g["rank"], g["gene"])],
                         columns=["cancer", "rank", "gene"])
    elif kind == "resc":
        t = pd.read_csv(RESC / fname)
        if "rank" not in t:
            t["rank"] = t.groupby("cancer").cumcount() + 1
    else:
        t = pd.read_csv(DERIVED / fname)
    t = t[t["rank"] <= 10][["cancer", "rank", "gene"]].copy()
    t["gene"] = t["gene"].astype(str).str.strip().str.upper()
    return t


def boot_ci(vec: np.ndarray):
    rng = np.random.default_rng(BOOT_SEED)
    n = len(vec)
    idx = rng.integers(0, n, size=(BOOT_B, n))
    m = vec[idx].mean(axis=1)
    return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def main():
    refs = load_references()
    for name, (g, e, nrows) in refs.items():
        print(f"reference {name}: {nrows} catalogue rows, {len(g)} genes, {len(e)} gene-cancer pairs")
    rows, vec = [], {}
    for sysname, (kind, fname) in SYSTEMS.items():
        t0 = load_system(kind, fname)
        for refname, (any_genes, exact, _) in refs.items():
            t = t0.copy()
            t["any"] = t["gene"].isin(any_genes)
            t["exact"] = [(g, c) in exact for g, c in zip(t["gene"], t["cancer"])]
            for metric, label in (("any", "any_cancer_P@10"), ("exact", "exact_cancer_concordance")):
                per = t.groupby("cancer")[metric].mean().reindex(CANCERS)
                v = per.to_numpy()
                lo, hi = boot_ci(v)
                rows.append(dict(system=sysname, reference=refname, metric=label,
                                 hits=int(t[metric].sum()), n=len(t),
                                 percent=round(100 * t[metric].sum() / len(t), 1),
                                 mean_P_over_cancers=float(v.mean()), boot_ci_lo=lo, boot_ci_hi=hi))
                vec[(sysname, refname, label)] = v
    summary = pd.DataFrame(rows)
    paired = []
    for (sysname, refname, label), v in vec.items():
        if sysname == "OncoDigger":
            continue
        d = vec[("OncoDigger", refname, label)] - v
        lo, hi = boot_ci(d)
        paired.append(dict(system=sysname, reference=refname, metric=label,
                           mean_diff_pp=round(100 * float(d.mean()), 1),
                           boot_ci_lo_pp=round(100 * lo, 1), boot_ci_hi_pp=round(100 * hi, 1),
                           cancers_oncodigger_better=int((d > 0).sum()), ties=int((d == 0).sum()),
                           cancers_oncodigger_worse=int((d < 0).sum())))
    paired = pd.DataFrame(paired)
    od = load_system("canon", None)
    pc = pd.DataFrame({"cancer": CANCERS})
    for refname, (any_genes, exact, _) in refs.items():
        a = od.assign(a=od["gene"].isin(any_genes), e=[(g, c) in exact for g, c in zip(od["gene"], od["cancer"])])
        agg = a.groupby("cancer")[["a", "e"]].sum().astype(int)
        pc[f"any_hits_{refname}"] = pc["cancer"].map(agg["a"])
        pc[f"exact_hits_{refname}"] = pc["cancer"].map(agg["e"])
    pc.to_csv(RESC / "statistical_support_reference_by_cancer_oncodigger.csv", index=False)
    summary.to_csv(RESC / "statistical_support_reference_summary.csv", index=False)
    paired.to_csv(RESC / "statistical_support_reference_paired_vs_oncodigger.csv", index=False)
    pd.set_option("display.width", 200)
    for label in ("any_cancer_P@10", "exact_cancer_concordance"):
        print(f"\n{label} (hits of 190)")
        print(summary[summary.metric == label].pivot(index="system", columns="reference", values="hits").to_string())
    print()
    print(paired[paired.reference == "statistical"].to_string(index=False))
    print(f"\nwrote {RESC / 'statistical_support_reference_summary.csv'}")


if __name__ == "__main__":
    main()
