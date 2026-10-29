"""
rank_aware_metrics.py -- nDCG, mean average precision and recall of the exact cancer-type IntOGen catalogue pairs, per cancer and averaged over the
19 cancers with a cancer-level bootstrap (19 cancers resampled, 2,000 draws, seed 42), for OncoDigger and the controls whose rankings are stored.

Relevance: a ranked (gene, cancer) pair is relevant if IntOGen lists that gene as a driver of that study cancer (the exact-concordance definition),
R_c is the number of such catalogue pairs for cancer c (1,740 in total).
  nDCG@k = DCG@k / IDCG@k with binary gains and a log2(i+1) discount, IDCG@k assuming min(k, R_c) relevant genes at the top.
  AP@k   = (1 / min(R_c, k)) * sum over ranks i <= k of precision@i * relevant_i.
  Recall@k = relevant genes in the top k / R_c.
The BM25-matched counting baseline stores only its top 10, so it is reported at k = 10; the others are reported at k = 10 and 50.
The pooled cancer-blind control has one list for all cancers.

Run    : python -X utf8 scripts/analysis/rank_aware_metrics.py
Output : data/derived/rank_aware_metrics.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "analysis"))
import interaction_all19 as ia  # noqa: E402

G, E = ia.load_intogen()
CANCERS = sorted({c for _, c in E})
R = {c: sum(1 for (_, cc) in E if cc == c) for c in CANCERS}


def systems():
    can = pd.read_csv(ROOT / "data" / "canonical" / "precision_at_k_detail.csv")
    pooled = pd.read_csv(ROOT / "data" / "derived" / "pooled_control_mutation.csv").sort_values("rank")
    pl = pooled.gene.tolist()
    out = {
        "OncoDigger (full formula)": {c: g.sort_values("rank").gene.tolist() for c, g in can.groupby("cancer")},
        "Pooled cancer-blind control": {c: pl for c in CANCERS},
        "Frequency-only baseline (aliases, full corpus)": {c: g.sort_values("rank").gene.tolist() for c, g in
                                                           pd.read_csv(ROOT / "data" / "derived" / "freq_baseline_top50_aliases.csv").groupby("cancer")},
        "BM25-matched flat counting (top 10 only)": {c: g.sort_values("rank").gene.tolist() for c, g in
                                                    pd.read_csv(ROOT / "data" / "derived" / "intogen_rescoring" / "matched_frequency_baseline_detail.csv").groupby("cancer")},
    }
    return out


def metrics(genes, cancer, k):
    rel = [1 if (g.upper(), cancer) in E else 0 for g in genes[:k]]
    rc = R[cancer]
    dcg = sum(r / np.log2(i + 2) for i, r in enumerate(rel))
    idcg = sum(1 / np.log2(i + 2) for i in range(min(k, rc)))
    ap, hits = 0.0, 0
    for i, r in enumerate(rel, 1):
        if r:
            hits += 1
            ap += hits / i
    return dcg / idcg, ap / min(rc, k), sum(rel) / rc, sum(rel)


def main() -> None:
    rng = np.random.default_rng(42)
    idx = rng.integers(0, len(CANCERS), size=(2000, len(CANCERS)))
    rows = []
    for name, lists in systems().items():
        for k in (10, 50):
            if "top 10 only" in name and k == 50:
                continue
            per = np.array([metrics(lists[c], c, k) for c in CANCERS])
            row = {"system": name, "k": k}
            for j, m in enumerate(("nDCG", "MAP", "recall_macro")):
                mean = per[:, j].mean()
                boots = per[idx, j].mean(axis=1)
                row[m] = round(100 * mean, 1)
                row[m + "_ci"] = f"{100 * np.percentile(boots, 2.5):.1f}-{100 * np.percentile(boots, 97.5):.1f}"
            row["recall_pooled"] = round(100 * per[:, 3].sum() / sum(R.values()), 1)
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "data" / "derived" / "rank_aware_metrics.csv", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
