"""Score counterfactual rankings: rank, support, relevance, top-10 membership, and the benchmark's two IntOGen metrics (published reference)."""
import os
import sys, json
from pathlib import Path
import pandas as pd
REPO=Path(os.environ.get("ONCODIGGER_PAPER_REPO") or Path(__file__).resolve().parents[3]); sys.path.insert(0,str(REPO/"scripts/analysis"))
import intogen_statistical_support_rescoring as base
anyg,exact,_=base.load_references()["published"]
canon=pd.read_csv(REPO/"data/canonical/precision_at_k_detail.csv")[["cancer","rank","gene"]]; canon=canon[canon["rank"]<=10].copy()
LAB={"breast":"Breast Cancer","lung":"Lung Cancer"}
def bench(cancer,top10):
    t=canon[canon.cancer!=LAB[cancer]]; add=pd.DataFrame({"cancer":LAB[cancer],"rank":range(1,11),"gene":top10}); t=pd.concat([t,add])
    t["gene"]=t.gene.str.upper(); a=t.gene.isin(anyg); e=[(g,c) in exact for g,c in zip(t.gene,t.cancer)]
    c=t[t.cancer==LAB[cancer]]
    return int(a.sum()),int(sum(e)),int(c.gene.isin(anyg).sum()),int(sum((g,LAB[cancer]) in exact for g in c.gene))
for cancer in ("breast","lung"):
    files=sorted(Path(".").glob(f"cf_{cancer}_*.json"))+(sorted(Path(".").glob("cf_baseline.json"))+sorted(Path(".").glob("cf_erbb2*.json")) if cancer=="breast" else [])
    if not files: continue
    print("\n=====",LAB[cancer])
    committed=canon[canon.cancer==LAB[cancer]].gene.tolist()
    print("committed top10:",committed)
    for f in files:
        d=json.load(open(f)); top=[o["gene"] for o in d[:10]]
        h=bench(cancer,top)
        print(f"\n{f.name}: top10={top}\n   same as committed top-10 (order): {top==committed}; same set: {set(top)==set(committed)}")
        print(f"   benchmark any-cancer {h[0]}/190, exact {h[1]}/190 | this cancer any {h[2]}/10, exact {h[3]}/10")
        for o in d[:14]: print(f"   {o['rank']:>2} {o['gene']:7} rel={o['relevance']:.3f} support={o['support']:>4} B={o['background']}")
