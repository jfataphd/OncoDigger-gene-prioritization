"""
formula_variants.py -- how sensitive are the top-10 gene lists to the form of the relevance score? (exploratory, read-only on features_raw.csv)

Uses the per-gene components extracted by extract_features.py (same pools, vocabulary and options as the paper; the recomputed baseline order
reproduces data/derived/query_planes/mutation_vocab_corrected.csv exactly). Every variant re-ranks the SAME candidate genes of the SAME pool,
so only the scoring changes.

IMPORTANT. IntOGen agreement is reported as a sensitivity measure only. Choosing a formula because it scores higher against IntOGen would be
tuning on the benchmark, which the paper states was not done. Any variant adopted later needs evidence that does not use IntOGen.

Variants
  baseline            R = W * ln(1+E) * ln(1+S)
  drop one factor     W*ln(1+E) | W*ln(1+S) | ln(1+E)*ln(1+S)
  single factor       W | S (flat count) | ln(1+E)
  rank decay in W     none (=S), 1/r, 1/log2(r+1), 1/sqrt(r) [baseline], 1/r^2   (each inside the full formula)
  enrichment form     ln(1+E) [baseline], ln(E), -log10 hypergeometric (Fisher) p, G statistic, PMI-style S*ln(lift)
  IDF style           S*ln(C/B), W*ln(C/B), W*ln(1+E) with B floor
  de-duplicated       (W/S) * ln(1+E) * ln(1+S)   (mean rank weight; support enters once)
Notation: S supporting papers in the pool, N pool size, B background papers, C collection size, E = (S/N)/(B/C), r = 1-based pool rank.

Run: python -X utf8 formula_variants.py     Output: variant_summary.csv, variant_top10.csv
"""
import importlib.util
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
spec = importlib.util.spec_from_file_location("ir", REPO / "scripts" / "analysis" / "intogen_full_rescoring.py")
ir = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ir)
GA, PAIRS = ir.load_intogen()

f = pd.read_csv(HERE / "features_raw.csv")
f["Be"] = np.maximum(f.B, f.S)
f["E"] = (f.S / f.N) / (f.Be / f.C)
rk = [np.array([int(x) for x in s.split(";")], float) for s in f.ranks]


def wsum(fn):
    return np.array([fn(r).sum() for r in rk])


W = f.W.values
S = f.S.values.astype(float)
E = f.E.values
lnE1 = np.log1p(E)
lnS1 = np.log1p(S)
Wd = {"none": wsum(lambda r: np.ones_like(r)), "1/r": wsum(lambda r: 1 / r), "1/log2(r+1)": wsum(lambda r: 1 / np.log2(r + 1)),
      "1/sqrt(r)": W, "1/r^2": wsum(lambda r: 1 / r ** 2)}
fisher = -hypergeom.logsf(f.S.values - 1, f.C.values, f.Be.values, f.N.values) / np.log(10)
Gstat = 2 * (S * np.log(E) + 0) if False else 2 * (S * np.log(np.maximum(E, 1e-12)))
lift = np.log(np.maximum(E, 1e-12))

V = {
    "baseline W*ln(1+E)*ln(1+S)": W * lnE1 * lnS1,
    "drop ln(1+S):  W*ln(1+E)": W * lnE1,
    "drop enrichment:  W*ln(1+S)": W * lnS1,
    "drop W:  ln(1+E)*ln(1+S)": lnE1 * lnS1,
    "W only": W,
    "S only (flat count)": S,
    "ln(1+E) only": lnE1,
    "decay none (W=S): S*ln(1+E)*ln(1+S)": Wd["none"] * lnE1 * lnS1,
    "decay 1/r": Wd["1/r"] * lnE1 * lnS1,
    "decay 1/log2(r+1)": Wd["1/log2(r+1)"] * lnE1 * lnS1,
    "decay 1/r^2": Wd["1/r^2"] * lnE1 * lnS1,
    "enrichment ln(E) instead of ln(1+E)": W * np.where(lift > 0, lift, 0) * lnS1,
    "enrichment -log10 Fisher p (W*p)": W * fisher,
    "enrichment -log10 Fisher p alone": fisher,
    "enrichment G statistic alone": Gstat,
    "S*ln(lift) (PMI-style)": S * np.where(lift > 0, lift, 0),
    "IDF style S*ln(C/B)": S * np.log(f.C.values / f.Be.values),
    "IDF style W*ln(C/B)": W * np.log(f.C.values / f.Be.values),
    "de-duplicated (W/S)*ln(1+E)*ln(1+S)": (W / S) * lnE1 * lnS1,
}

rows, tops = [], []
base_top = None
for name, sc in V.items():
    f["_s"] = sc
    t = {}
    for c, g in f.groupby("cancer", sort=False):
        t[c] = g.assign(_o=g.baseline_rank).sort_values(["_s", "_o"], ascending=[False, True], kind="stable").gene.head(10).tolist()
    if base_top is None:
        base_top = t
    d10 = pd.DataFrame([{"cancer": c, "gene": x} for c, v in t.items() for x in v])
    sc_ = ir.score(d10, GA, PAIRS)
    per = sc_.groupby("cancer")[["intogen_any", "intogen_exact"]].sum().reindex(ir.STUDY_ORDER)
    jb = np.mean([len(set(t[c]) & set(base_top[c])) / len(set(t[c]) | set(base_top[c])) for c in t])
    jp = np.mean([len(set(t[a]) & set(t[b])) / len(set(t[a]) | set(t[b])) for a, b in combinations(t, 2)])
    nsh = np.mean([len(set(t[c]) & set(base_top[c])) for c in t])
    rows.append({"variant": name, "any_of_190": int(per.intogen_any.sum()), "exact_of_190": int(per.intogen_exact.sum()),
                 "mean_top10_overlap_with_baseline": round(nsh, 2), "mean_top10_jaccard_with_baseline": round(jb, 3),
                 "mean_pairwise_jaccard_between_cancers": round(jp, 3), "per_cancer_exact": per.intogen_exact.tolist()})
    for c, v in t.items():
        for i, x in enumerate(v, 1):
            tops.append({"variant": name, "cancer": c, "rank": i, "gene": x})
out = pd.DataFrame(rows)
# paired cancer-bootstrap difference in exact concordance versus baseline
rng = np.random.default_rng(42)
idx = rng.integers(0, 19, size=(2000, 19))
b = np.array(out.per_cancer_exact[0], float)
lo, hi = [], []
for pc in out.per_cancer_exact:
    d = (np.array(pc, float) - b) / 10 * 100
    bs = np.array([d[i].mean() for i in idx])
    lo.append(round(np.percentile(bs, 2.5), 1)); hi.append(round(np.percentile(bs, 97.5), 1))
out["exact_diff_vs_baseline_pts_CI"] = [f"{(e - out.exact_of_190[0]) / 1.9:+.1f} ({l:+.1f} to {h:+.1f})" for e, l, h in zip(out.exact_of_190, lo, hi)]
out = out.drop(columns="per_cancer_exact")
out.to_csv(HERE / "variant_summary.csv", index=False)
pd.DataFrame(tops).to_csv(HERE / "variant_top10.csv", index=False)
pd.set_option("display.width", 250)
print(out.to_string(index=False))
