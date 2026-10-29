"""
build_primary_tables.py — primary-ranking tables from an OncoDigger ranking file.

Input
  --ranking   CSV with columns cancer, rank, gene (ranks 1-50), e.g.
              data/derived/query_planes/mutation.csv (from run_query_planes.py)
  data/derived/pooled_control_mutation.csv              the pooled-corpus control ranking

Outputs (--out-dir, default data/canonical)
  precision_at_k_detail.csv            per cancer x rank: gene and rank tier (the canonical ranking read by every later script;
                                       the file name is historical, it holds no precision values)
  ../derived/corpus_specificity_clean.csv  overlap with the pooled-corpus global top-10

Run:  python -X utf8 scripts/analysis/build_primary_tables.py --ranking data/derived/query_planes/mutation.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
# global control: the same pipeline on all 19 corpora pooled (pooled_corpus_control.py)
GLOBAL_CONTROL = ROOT / "data" / "derived" / "pooled_control_mutation.csv"


def tier_of(k: int) -> str:
    return "established" if k <= 10 else "emerging" if k <= 20 else "discovery" if k <= 30 else "extended"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ranking", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, default=ROOT / "data" / "canonical")
    args = ap.parse_args()

    rk = pd.read_csv(args.ranking).sort_values(["cancer", "rank"])
    g = pd.read_csv(GLOBAL_CONTROL).sort_values("rank")
    GLOBAL_TOP10 = g.loc[g["rank"] <= 10, "gene"].tolist()

    # canonical ranking with rank tiers
    rows = []
    for c, g in rk.groupby("cancer", sort=False):
        for k, gene in enumerate(g["gene"], 1):
            rows.append({"cancer": c, "rank": k, "gene": gene, "tier": tier_of(k)})
    detail = pd.DataFrame(rows)

    top10 = {c: dict(zip(g["gene"], g["rank"])) for c, g in rk[rk["rank"] <= 10].groupby("cancer")}

    spec = []
    for c, genes in top10.items():
        s = set(genes)
        specific = sorted(s - set(GLOBAL_TOP10))
        spec.append({"cancer": c, "overlap_with_global": len(s & set(GLOBAL_TOP10)),
                     "cancer_specific_in_top10": len(specific),
                     "cancer_specific_genes": "; ".join(specific)})
    spec = pd.DataFrame(spec).sort_values("cancer")

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    detail.to_csv(out / "precision_at_k_detail.csv", index=False)
    derived = ROOT / "data" / "derived" if out == ROOT / "data" / "canonical" else out
    spec.to_csv(derived / "corpus_specificity_clean.csv", index=False)
    print(f"{len(detail)} ranked pairs; mean cancer-specific genes in the top 10 {spec['cancer_specific_in_top10'].mean():.1f}; "
          f"mean overlap with the pooled control {spec['overlap_with_global'].mean():.1f}")


if __name__ == "__main__":
    main()
