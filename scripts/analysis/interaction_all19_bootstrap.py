"""
interaction_all19_bootstrap.py -- pooled retrieval x vocabulary interaction and its case-bootstrap
95% CI over all 19 cancers (Results 5.8, Table 2, Abstract, Supplementary Figure S9).

The point estimate and CI quoted for the 19-cancer interaction previously had no generating script in
this repository (interaction_bootstrap_5cancer.py covers only its 5-cancer diagnostic sample). This
script applies the same case-bootstrap used there, unchanged (resample cancers with replacement,
2,000 replicates, seed 42, percentile bounds at sorted positions int(0.025*B) and int(0.975*B)), to
the per-cancer interaction terms in interaction_all19_final.csv.

Input:  data/derived/intogen_rescoring/interaction_all19_final.csv   (interaction_all19.py)
Output: data/derived/intogen_rescoring/interaction_all19_bootstrap_ci.csv

Run: python -X utf8 scripts/analysis/interaction_all19_bootstrap.py [path/to/interaction_all19_final.csv]
"""
from __future__ import annotations

import csv
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DER = ROOT / "data" / "derived" / "intogen_rescoring"


def bootstrap(rows: list[dict], prefix: str, n_boot: int = 2000, seed: int = 42) -> tuple[float, float, float]:
    cells = {k: {r["cancer"]: int(r[f"{k}_{prefix}"]) for r in rows} for k in ("cell1", "cell2", "cell3", "cell4")}
    cancers = [r["cancer"] for r in rows]
    rng = random.Random(seed)
    stats = []
    for _ in range(n_boot):
        sample = [rng.choice(cancers) for _ in cancers]
        c = {k: sum(cells[k][x] for x in sample) for k in cells}
        stats.append(100 * (c["cell4"] - c["cell2"] - c["cell3"] + c["cell1"]) / (len(sample) * 10))
    stats.sort()
    point = 100 * sum(cells["cell4"][x] - cells["cell2"][x] - cells["cell3"][x] + cells["cell1"][x]
                      for x in cancers) / (len(cancers) * 10)
    return point, stats[int(0.025 * n_boot)], stats[int(0.975 * n_boot)]


def main() -> None:
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else DER / "interaction_all19_final.csv"
    rows = list(csv.DictReader(open(src, encoding="utf-8")))
    assert len(rows) == 19, len(rows)
    out = []
    for metric, prefix in (("any_cancer_P@10", "any"), ("exact_concordance", "exact")):
        point, lo, hi = bootstrap(rows, prefix)
        per = [int(r[f"{prefix}_interaction"]) for r in rows]
        out.append({"metric": metric, "interaction_pp": round(point, 4), "ci_lo_pp": round(lo, 4),
                    "ci_hi_pp": round(hi, 4), "n_cancers_positive": sum(1 for v in per if v > 0),
                    "per_cancer_min": min(per), "per_cancer_max": max(per)})
        print(f"{metric}: interaction {point:+.1f} pp  95% case-bootstrap CI {lo:+.1f} to {hi:+.1f}  "
              f"positive in {sum(1 for v in per if v > 0)}/19 cancers (range {min(per):+d} to {max(per):+d})")
    if len(sys.argv) == 1:
        with open(DER / "interaction_all19_bootstrap_ci.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(out[0]))
            w.writeheader()
            w.writerows(out)
        print("Wrote", DER / "interaction_all19_bootstrap_ci.csv")


if __name__ == "__main__":
    main()
