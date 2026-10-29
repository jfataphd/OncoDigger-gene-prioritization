"""
gene_exclusion_readmission_test.py -- post hoc test of the three inherited full gene exclusions (FBN1, SGCG, WDHD1;
Supplementary Methods S1). Re-runs the primary `mutation mutations` ranking (same engine call as run_query_planes.py)
for all 19 cancers with the exclusion lifted, under three vocabularies:

  A  every term admitted (symbol, approved name, all aliases)
  B  as A with the bare symbol suppressed (control: identical to A, so the effect is not carried by the symbol)
  C  approved symbol and approved full name only (every alias of the three genes dropped)

and compares the top-50 lists with the committed run on the same corpora (--baseline).

Env: ONCODIGGER_ROOT = verified 2026-09-27 corpus root; PYTHONHASHSEED=0 (set by run_query_planes' own re-exec rule).
Usage: PYTHONHASHSEED=0 python -X utf8 scripts/analysis/gene_exclusion_readmission_test.py run A|B|C
       python -X utf8 scripts/analysis/gene_exclusion_readmission_test.py compare
Outputs: data/derived/gene_exclusion_test/mutation_<variant>.csv and summary.csv
"""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "derived" / "gene_exclusion_test"
THREE = {"FBN1", "SGCG", "WDHD1"}
KEEP_PHRASES = {"fibrillin 1", "sarcoglycan gamma", "wd repeat and hmg box dna binding protein 1"}


def load_module():
    spec = importlib.util.spec_from_file_location("rqp", ROOT / "scripts" / "analysis" / "run_query_planes.py")
    m = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(ROOT / "scripts" / "analysis"))
    spec.loader.exec_module(m)
    return m


def run(variant: str) -> None:
    m = load_module()
    m.GENE_EXCLUDE.difference_update(THREE)
    if variant == "B":
        m.SUPPRESS_SYMBOLS.update(THREE)
    if variant == "C":
        orig = m.suppress_ambiguous_symbols

        def patched(lex):
            out = orig(lex)
            single = {}
            for t, es in out["single"].items():
                kept = [e for e in es if not (e["id"] in THREE and t != e["id"].lower())]
                if kept:
                    single[t] = kept
            phr = {}
            for a, lst in out["phrases"].items():
                keep = [(term, e) for term, e in lst if not (e["id"] in THREE and term.replace("-", " ") not in KEEP_PHRASES)]
                if keep:
                    phr[a] = keep
            return {"single": single, "phrases": phr}

        m.suppress_ambiguous_symbols = patched
    OUT.mkdir(parents=True, exist_ok=True)
    tmp = OUT / f"_tmp_{variant}"
    sys.argv = [str(ROOT / "scripts" / "analysis" / "run_query_planes.py"), "--planes", "mutation",
                "--oncodigger-root", os.environ["ONCODIGGER_ROOT"], "--out-dir", str(tmp)]
    m.main()
    (tmp / "mutation.csv").replace(OUT / f"mutation_{variant}.csv")
    tmp.rmdir()


def compare() -> None:
    base = pd.read_csv(ROOT / "data" / "derived" / "query_planes" / "mutation.csv")
    canon = pd.read_csv(ROOT / "data" / "canonical" / "precision_at_k_detail.csv")
    rows = []
    for v in "ABC":
        d = pd.read_csv(OUT / f"mutation_{v}.csv")
        mm = base.merge(d, on=["cancer", "rank"], suffixes=("_b", "_v"))
        r = {"variant": v}
        for k in (10, 50):
            s = mm[mm["rank"] <= k]
            r[f"top{k}_slots_changed"] = int((s.gene_b != s.gene_v).sum())
        t50 = d[(d["rank"] <= 50) & d.gene.isin(THREE)]
        r["three_gene_pairs_in_top50"] = len(t50)
        r["three_gene_pairs_in_top10"] = int((t50["rank"] <= 10).sum())
        r["cancers_identical_top50_set"] = int(sum(
            set(base[(base.cancer == c) & (base["rank"] <= 50)].gene) == set(d[(d.cancer == c) & (d["rank"] <= 50)].gene)
            for c in base.cancer.unique()))
        rows.append(r)
    pd.DataFrame(rows).to_csv(OUT / "summary.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    run(sys.argv[2]) if sys.argv[1] == "run" else compare()
