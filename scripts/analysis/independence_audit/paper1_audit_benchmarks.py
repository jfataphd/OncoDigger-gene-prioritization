"""Recompute Paper 1 benchmark numbers under stricter, knowledge-independent-as-possible IntOGen references (sensitivity analysis only).

References (gene-level 'any cancer' G*, and gene-cancer 'exact' E*), all from the released 2024.09.20 unfiltered table:
  0 published        : catalog rows (FILTER==PASS)                                  [reproduces 177/190 and 144/190]
  1 stat_in_catalog  : catalog rows with q<=0.05 and >=2 significant methods         (drops CGC-waiver and CGC-rescue rows)
  2 classA           : catalog rows of audit class A (statistic + no CGC-only filter exemption)
  3 stat_any_filter  : ANY unfiltered row with q<=0.05 and >=2 significant methods   (statistical call before IntOGen's vetting filters)
  4 no_lit_gate      : catalog rows plus the rows removed only by the 'Lack of literature evidence' filter (counterfactual: CancerMine gate removed)
Not an exact IntOGen rerun without CGC/CancerMine: the pipeline was not re-executed. Method weights (CGC-calibrated) cannot be removed.
CI: case bootstrap over the 19 cancers (2000 resamples, seed 42), as in Paper 1.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paper1_audit_main as M  # noqa: E402

P1, OUT = M.P1, M.OUT
SYSTEMS = {
    "OncoDigger": ("canon", None),
    "pooled_control_cancer_blind": ("pooled", None),
    "CancerMine": ("file", "cancermine_top10.csv"),
    "ChatGPT_free_tier": ("file", "chatgpt_top10.csv"),
    "DISEASES_2.0": ("file", "diseases_top10.csv"),
    "PubTator_3.0_v2": ("file", "pubtator_v2_top10.csv"),
    "frequency_baseline_aliases": ("freq", "freq_baseline_top50_aliases.csv"),
    "frequency_baseline_primary": ("freq", "freq_baseline_top50_primary.csv"),
    "matched_frequency_only": ("resc", "matched_frequency_baseline_detail.csv"),
    "ablation_full_pipeline": ("resc", "ablation_full_pipeline_detail.csv"),
    "ablation_bm25_vetted_no_longest_match_fix": ("resc", "ablation_bm25_vetted_no_longest_match_fix_detail.csv"),
    "ablation_bm25_raw_uncleaned_vocab": ("resc", "ablation_bm25_raw_uncleaned_vocab_detail.csv"),
    "ablation_frequency_only_vetted_aliases": ("resc", "ablation_frequency_only_vetted_aliases_detail.csv"),
    "ablation_frequency_only_primary_symbols": ("resc", "ablation_frequency_only_primary_symbols_detail.csv"),
}


def systems_tables(cancers):
    out = {}
    d = P1 / "data/derived"
    for name, (kind, f) in SYSTEMS.items():
        if kind == "canon":
            t = pd.read_csv(P1 / "data/canonical/precision_at_k_detail.csv")
            t = t[t["rank"] <= 10][["cancer", "rank", "gene"]]
        elif kind == "pooled":
            g = pd.read_csv(d / "pooled_control_mutation.csv").head(10)
            t = pd.DataFrame([(c, r, x) for c in cancers for r, x in zip(g["rank"], g.gene)], columns=["cancer", "rank", "gene"])
        elif kind == "resc":
            t = pd.read_csv(d / "intogen_rescoring" / f)
            if "rank" not in t:
                t["rank"] = t.groupby("cancer").cumcount() + 1
            t = t[t["rank"] <= 10][["cancer", "rank", "gene"]]
        else:
            t = pd.read_csv(d / f)
            t = t[t["rank"] <= 10][["cancer", "rank", "gene"]]
        t = t.copy()
        t["gene"] = t.gene.str.upper()
        out[name] = t
    return out


def boot(vec, B=2000, seed=42):
    rng = np.random.default_rng(seed)
    n = len(vec)
    idx = rng.integers(0, n, size=(B, n))
    m = vec[idx].mean(axis=1)
    return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def main():
    P, C, D, comp, unf, CMAP = M.main()
    cancers = list(CMAP)
    p = unf.is_pass
    lack = unf.FILTER == "Lack of literature evidence"
    refs = {
        "0_published": unf[p],
        "1_stat_in_catalog": unf[p & unf.statistical_call_noncgc_rule],
        "2_classA": unf[unf.row_class == "A"],
        "3_stat_any_filter": unf[unf.statistical_call_noncgc_rule],
        "4_no_literature_gate": unf[p | lack],
    }
    sets = {}
    for k, df in refs.items():
        sets[k] = (set(df.SYMBOL), set(zip(df.SYMBOL, df.study_cancer.dropna().reindex(df.index).fillna(""))) if False else set(zip(df.SYMBOL, df.study_cancer)))
    S = systems_tables(cancers)
    rows, percancer = [], []
    vec = {}
    for sysname, t in S.items():
        for rk, (gset, eset) in sets.items():
            t = t.copy()
            t["any"] = t.gene.isin(gset)
            t["exact"] = [(g, c) in eset for g, c in zip(t.gene, t.cancer)]
            for metric in ("any", "exact"):
                per = t.groupby("cancer")[metric].mean().reindex(cancers)
                v = per.to_numpy()
                lo, hi = boot(v)
                rows.append(dict(system=sysname, reference=rk, metric="any_cancer_P@10" if metric == "any" else "exact_cancer_concordance", hits=int(t[metric].sum()), n=len(t),
                                 mean_P_over_cancers=float(v.mean()), boot_ci_lo=lo, boot_ci_hi=hi, cancers_at_100pct=int((v == 1).sum()), cancers_below_50pct=int((v < 0.5).sum())))
                vec[(sysname, rk, metric)] = v
                if sysname == "OncoDigger":
                    for c, x in per.items():
                        percancer.append(dict(cancer=c, reference=rk, metric=metric, precision=float(x), hits=int(round(x * 10))))
    R = pd.DataFrame(rows)
    # paired differences vs OncoDigger
    pr = []
    for (sysname, rk, metric), v in vec.items():
        if sysname == "OncoDigger":
            continue
        d = vec[("OncoDigger", rk, metric)] - v
        lo, hi = boot(d)
        pr.append(dict(system=sysname, reference=rk, metric=metric, mean_diff_oncodigger_minus_system=float(d.mean()), boot_ci_lo=lo, boot_ci_hi=hi,
                       cancers_oncodigger_better=int((d > 0).sum()), ties=int((d == 0).sum()), cancers_oncodigger_worse=int((d < 0).sum())))
    R.to_csv(OUT / "benchmark_recalculation.csv", index=False)
    pd.DataFrame(pr).to_csv(OUT / "benchmark_recalculation_paired_vs_oncodigger.csv", index=False)
    pd.DataFrame(percancer).pivot_table(index="cancer", columns=["metric", "reference"], values="precision").to_csv(OUT / "benchmark_recalculation_by_cancer_oncodigger.csv")
    od = R[R.system == "OncoDigger"]
    pd.set_option("display.width", 220)
    print(od[["reference", "metric", "hits", "n", "mean_P_over_cancers", "boot_ci_lo", "boot_ci_hi", "cancers_below_50pct"]].to_string(index=False))
    sub = R[R.metric == "exact_cancer_concordance"].pivot(index="system", columns="reference", values="hits")
    print(sub.to_string())
    sub = R[R.metric == "any_cancer_P@10"].pivot(index="system", columns="reference", values="hits")
    print(sub.to_string())


if __name__ == "__main__":
    main()
