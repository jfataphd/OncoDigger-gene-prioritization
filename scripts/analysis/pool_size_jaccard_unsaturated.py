"""
pool_size_jaccard_unsaturated.py -- pool-size Jaccard restricted to cancers whose pool is not exhausted.

For Cervical Cancer and Myeloma only 861 and 817 abstracts contain a query word, so every pool size above that returns the same papers and
the same top 10. Those two cancers contribute Jaccard = 1 for the comparisons that involve pools of 1,000 or more, which raises the all-19
means. This script recomputes the means over the other 17 cancers from the stored per-size lists (no engine rerun) and states, for each
cancer, whether the pools of 1,000, 1,500 and 2,000 gave identical top-10 lists.

Input : data/derived/pool_size_top10.csv (pool_size_jaccard.py)
Output: data/derived/pool_size_jaccard_17cancers.csv
Run   : python -X utf8 scripts/analysis/pool_size_jaccard_unsaturated.py
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
d = pd.read_csv(ROOT / "data" / "derived" / "pool_size_top10.csv")
top = {(r.pool_size, r.cancer): r.top10.split("|") for r in d.itertuples()}
cancers = sorted(d.cancer.unique())
SATURATED = ["Cervical Cancer", "Myeloma"]


def jac(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if a | b else 1.0


def means(cs):
    pairs = [(250, 1000), (500, 1000), (750, 1000), (1500, 1000), (2000, 1000)]
    adj = [(250, 500), (500, 750), (750, 1000), (1000, 1500), (1500, 2000)]
    out = {f"{a} vs {b}": sum(jac(top[a, c], top[b, c]) for c in cs) / len(cs) for a, b in pairs}
    for a, b in adj:
        out[f"{a}/{b} adjacent"] = sum(jac(top[a, c], top[b, c]) for c in cs) / len(cs)
    out["mean of adjacent"] = sum(v for k, v in out.items() if k.endswith("adjacent")) / 5
    return out


rows = []
for label, cs in (("all 19", cancers), ("17 (Cervical Cancer and Myeloma excluded)", [c for c in cancers if c not in SATURATED])):
    for k, v in means(cs).items():
        rows.append({"cancers": label, "comparison": k, "mean_top10_jaccard": round(v, 3)})
out = pd.DataFrame(rows)
out.to_csv(ROOT / "data" / "derived" / "pool_size_jaccard_17cancers.csv", index=False)
print(out.pivot(index="comparison", columns="cancers", values="mean_top10_jaccard").to_string())
print()
for c in cancers:
    same = top[1000, c] == top[1500, c] == top[2000, c]
    if same:
        print("identical top 10 at 1,000, 1,500 and 2,000:", c)
