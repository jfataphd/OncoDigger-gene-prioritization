"""
driver_mutation_sensitivity.py -- wording sensitivity of the primary ranking (Supplementary Results S3): the secondary query
`driver mutation driver mutations` against the primary query `mutation mutations`, for all 19 cancers.

Reads data/derived/query_planes/driver_mutation.csv (run_query_planes.py --planes driver_mutation) and the canonical primary ranking
data/canonical/precision_at_k_detail.csv, scores both against IntOGen (any-cancer and exact cancer-type concordance, the same reference,
mapping and 2,000-draw seed-42 cancer bootstrap as intogen_full_rescoring.py) and reports the mean top-10 Jaccard similarity of the two
lists per cancer (with range) and the paired t-test of each metric.
Output: data/derived/intogen_rescoring/driver_mutation_sensitivity.csv (one row per metric) and a printed summary.
Optional arguments: <driver_mutation.csv> <primary.csv> to score other files (used to check the earlier committed run).
Run: python -X utf8 scripts/analysis/driver_mutation_sensitivity.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "analysis"))
import intogen_full_rescoring as ir  # noqa: E402

DRV = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "data" / "derived" / "query_planes" / "driver_mutation.csv"
PRI = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "data" / "canonical" / "precision_at_k_detail.csv"
OUT = ROOT / "data" / "derived" / "intogen_rescoring" / "driver_mutation_sensitivity.csv"


def top10(path):
    d = pd.read_csv(path)
    d = d[d["rank"] <= 10]
    return {c: list(g.sort_values("rank").gene) for c, g in d.groupby("cancer")}


def main() -> None:
    ga, pairs = ir.load_intogen()
    a, b = top10(PRI), top10(DRV)
    cancers = sorted(a)
    assert sorted(b) == cancers and all(len(a[c]) == 10 and len(b[c]) == 10 for c in cancers)
    hit_any = lambda L, c: sum(g in ga for g in L[c]) / 10
    hit_ex = lambda L, c: sum((g, c) in pairs for g in L[c]) / 10
    jac = np.array([len(set(a[c]) & set(b[c])) / len(set(a[c]) | set(b[c])) for c in cancers])
    rows = []
    for name, f in (("any_cancer_P@10", hit_any), ("exact_cancer_concordance", hit_ex)):
        p = np.array([f(a, c) for c in cancers]); q = np.array([f(b, c) for c in cancers])
        pc, qc = (p * 10).round().astype(int), (q * 10).round().astype(int)    # integer hit counts, so that ties are exact
        _m, lo, hi, _ = ir.paired_bootstrap_ci(q, p)
        rows.append({"metric": name, "primary_hits": int(pc.sum()), "driver_query_hits": int(qc.sum()),
                     "primary_pct": 100 * p.mean(), "driver_query_pct": 100 * q.mean(), "diff_pp": 100 * (q - p).mean(),
                     "ci_lo_pp": 100 * lo, "ci_hi_pp": 100 * hi, "p_ttest": stats.ttest_rel(qc, pc).pvalue,
                     "p_wilcoxon": stats.wilcoxon(qc, pc).pvalue if (qc != pc).any() else float("nan"),
                     "mean_top10_jaccard": jac.mean(), "jaccard_min": jac.min(), "jaccard_max": jac.max()})
    res = pd.DataFrame(rows)
    if DRV == ROOT / "data" / "derived" / "query_planes" / "driver_mutation.csv":
        res.to_csv(OUT, index=False)
    pd.set_option("display.width", 220, "display.max_columns", 20)
    print(res.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
