"""
intogen_precision_bootstrap_ci.py -- manuscript1.1 (IntOGen-primary draft) ONLY. Case-bootstrap
95% CI of Precision@k (19 cancers resampled with replacement, 2000 replicates, seed 42, for
each rank k=1-50), scored against IntOGen's any-cancer driver set.

Input:  data/derived/intogen_rescoring/precision_at_k_detail_intogen.csv
Output: data/derived/intogen_rescoring/precision_at_k_bootstrap_ci_intogen.csv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DER = ROOT / "data" / "derived" / "intogen_rescoring"
N_BOOT, SEED = 2000, 42
STUDY_ORDER = ['Breast Cancer', 'Lung Cancer', 'Colorectal Cancer', 'Prostate Cancer', 'Melanoma', 'Bladder Cancer',
               'Kidney Cancer', 'Pancreatic Cancer', 'Liver Cancer', 'Stomach Cancer', 'Esophageal Cancer',
               'Ovarian Cancer', 'Endometrial Cancer', 'Cervical Cancer', 'Thyroid Cancer', 'Brain Cancer',
               'Leukemia', 'Lymphoma', 'Myeloma']


def main():
    d = pd.read_csv(DER / "precision_at_k_detail_intogen.csv")
    # Canonical bootstrap shared with intogen_statistical_support_rescoring.boot_ci: study cancer order, fresh
    # seed for each statistic, one vectorised draw of (N_BOOT, 19) indices, so k=10 equals the Table 3 interval.
    cancers = STUDY_ORDER
    n = len(cancers)

    rows = []
    for k in range(1, 51):
        sub = d[d["rank"] == k].set_index("cancer").loc[cancers, "precision_at_k"].values
        idx = np.random.default_rng(SEED).integers(0, n, size=(N_BOOT, n))
        boots = sub[idx].mean(axis=1)
        lo, hi = np.percentile(boots, [2.5, 97.5])
        rows.append({"k": k, "mean_p_at_k": sub.mean(), "ci_low_2.5": lo, "ci_high_97.5": hi,
                     "n_cancers": n, "n_bootstrap": N_BOOT})

    out = DER / "precision_at_k_bootstrap_ci_intogen.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"Wrote {out}")
    for k in (5, 10, 15, 20, 30, 40):
        r = [r for r in rows if r["k"] == k][0]
        print(f"  k={k}: {r['mean_p_at_k']:.1%} (95% CI {r['ci_low_2.5']:.1%}-{r['ci_high_97.5']:.1%})")


if __name__ == "__main__":
    main()
