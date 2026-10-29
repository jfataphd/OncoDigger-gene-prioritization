"""
proximity_vs_matched_counting.py -- paired comparison of the proximity-adjusted ranking (w = 3) with the matched frequency-only control of Results 5.8
(BM25-retrieved papers and cleaned vocabulary held fixed, genes ranked by the number of supporting papers).

Inputs : data/derived/proximity/ranking_w3.csv and ranking_w0.csv (proximity_boost_run.py), data/derived/intogen_rescoring/matched_frequency_baseline_detail.csv
Output : data/derived/proximity/vs_matched_counting.csv  (exact and any-cancer top-10 counts, paired difference in points, 95% cancer-bootstrap CI,
         cancers better / equal / worse). Same IntOGen reference, cancer mapping and bootstrap (2,000 draws, seed 42) as proximity_summary.py.
Run    : python -X utf8 scripts/analysis/proximity_vs_matched_counting.py
"""
import importlib.util
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "data" / "derived" / "proximity"


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / "analysis" / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


ir = load("intogen_full_rescoring")
GA, PAIRS = ir.load_intogen()
ORDER = ir.STUDY_ORDER


def per_cancer(path):
    d = pd.read_csv(path)
    d["gene"] = d.gene.astype(str).str.upper()
    t = ir.score(d[d["rank"] <= 10][["cancer", "gene"]], GA, PAIRS)
    return t.groupby("cancer")[["intogen_any", "intogen_exact"]].sum().reindex(ORDER)


prox = per_cancer(D / "ranking_w3.csv")
frq = per_cancer(ROOT / "data" / "derived" / "intogen_rescoring" / "matched_frequency_baseline_detail.csv")
rows = []
for metric in ("intogen_exact", "intogen_any"):
    m, lo, hi, _ = ir.paired_bootstrap_ci(prox[metric].values / 10.0, frq[metric].values / 10.0)
    rows.append({"metric": metric, "proximity_w3": int(prox[metric].sum()), "matched_counting": int(frq[metric].sum()),
                 "diff_pp": round(100 * m, 1), "ci_lo": round(100 * lo, 1), "ci_hi": round(100 * hi, 1),
                 "cancers_better": int((prox[metric] > frq[metric]).sum()), "cancers_equal": int((prox[metric] == frq[metric]).sum()),
                 "cancers_worse": int((prox[metric] < frq[metric]).sum())})
out = pd.DataFrame(rows)
out.to_csv(D / "vs_matched_counting.csv", index=False)
print(out.to_string(index=False))
