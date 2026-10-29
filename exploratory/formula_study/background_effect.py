"""background_effect.py -- effect of the engine's background count on the top-10 lists (exploratory; reads features_raw.csv, background_union.csv).
Compares the baseline score with the same score using B_union (all active aliases) in place of the engine's B. Output: background_effect.csv, prints summary."""
import importlib.util
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent; REPO = HERE.parents[1]
spec = importlib.util.spec_from_file_location("ir", REPO / "scripts" / "analysis" / "intogen_full_rescoring.py")
ir = importlib.util.module_from_spec(spec); spec.loader.exec_module(ir); GA, PAIRS = ir.load_intogen()
f = pd.read_csv(HERE / "features_raw.csv").merge(pd.read_csv(HERE / "background_union.csv"), on=["cancer", "gene"])
f["Be"] = np.maximum(f.B, f.S); f["Bu"] = np.maximum(f.B_union, f.S)
def R(B): return f.W * np.log1p((f.S / f.N) / (B / f.C)) * np.log1p(f.S)
f["R0"], f["R1"] = R(f.Be), R(f.Bu)
print("genes where engine B < B_union:", int((f.Be < f.Bu).sum()), "of", len(f), "| median B_union/B:", round(float((f.Bu / f.Be).median()), 3),
      "| B_union/B > 2:", int((f.Bu / f.Be > 2).sum()))
t = {}
for col in ("R0", "R1"):
    t[col] = {c: g.sort_values([col, "baseline_rank"], ascending=[False, True]).gene.head(10).tolist() for c, g in f.groupby("cancer", sort=False)}
for col in t:
    d10 = pd.DataFrame([{"cancer": c, "gene": x} for c, v in t[col].items() for x in v]); s = ir.score(d10, GA, PAIRS)
    print(col, "any", int(s.intogen_any.sum()), "exact", int(s.intogen_exact.sum()))
sw = [(c, sorted(set(t["R0"][c]) - set(t["R1"][c])), sorted(set(t["R1"][c]) - set(t["R0"][c]))) for c in t["R0"] if set(t["R0"][c]) != set(t["R1"][c])]
print("cancers with a changed top 10:", len(sw), "| slots changed:", sum(len(a) for _, a, _ in sw))
for c, a, b in sw: print(" ", c, "left", a, "entered", b)
top50 = {col: {c: g.sort_values([col, "baseline_rank"], ascending=[False, True]).gene.head(50).tolist() for c, g in f.groupby("cancer", sort=False)} for col in ("R0", "R1")}
print("top-50 slots replaced:", sum(len(set(top50["R0"][c]) - set(top50["R1"][c])) for c in top50["R0"]), "of 950")
f[["cancer", "gene", "S", "B", "B_union", "basis"]].to_csv(HERE / "background_effect.csv", index=False)
d = f[f.basis == "phrase text scan"]; print("phrase-only genes:", len(d), "| in baseline top 50:", int((d.baseline_rank <= 50).sum()))
