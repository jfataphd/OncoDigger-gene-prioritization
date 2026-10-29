"""
apply_vocab_correction_to_ranking.py -- builds the vocabulary-corrected full ranking that feeds the
canonical tables (Methods 3.3).

The canonical gene rankings used for the headline tables are the output of the main pipeline
(data/derived/query_planes/mutation.csv, run_query_planes.py) with the MET (Lung), REST (Thyroid) and
FH (Kidney) false-positive correction applied for the three affected cancers
(vocab_corrected_rerun.py + vocab_corrected_rescoring.py -> vocab_corrected_rankings.csv). All other
cancers are taken unchanged from the main pipeline. Previously this substitution was made by hand; this
script makes it a recorded, repeatable step.

Inputs:
  data/derived/query_planes/mutation.csv
  data/derived/intogen_rescoring/vocab_corrected_rankings.csv
Output:
  data/derived/query_planes/mutation_vocab_corrected.csv   (cancer, rank, gene; 19 cancers x 50 ranks)

Prints every position at which the corrected ranking differs from the uncorrected one.

Run: python -X utf8 scripts/analysis/apply_vocab_correction_to_ranking.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
MUT = ROOT / "data" / "derived" / "query_planes" / "mutation.csv"
CORR = ROOT / "data" / "derived" / "intogen_rescoring" / "vocab_corrected_rankings.csv"
OUT = ROOT / "data" / "derived" / "query_planes" / "mutation_vocab_corrected.csv"


def main() -> None:
    base = pd.read_csv(MUT)
    corr = pd.read_csv(CORR)
    assert set(corr["cancer"]) <= set(base["cancer"]), set(corr["cancer"]) - set(base["cancer"])
    assert base["cancer"].nunique() == 19 and len(base) == 950
    merged = pd.concat([base[~base["cancer"].isin(corr["cancer"])], corr[["cancer", "rank", "gene"]]])
    order = {c: i for i, c in enumerate(base["cancer"].drop_duplicates())}
    merged = merged.assign(_o=merged["cancer"].map(order)).sort_values(["_o", "rank"]).drop(columns="_o")
    assert len(merged) == 950 and merged.groupby("cancer")["rank"].apply(lambda r: list(r) == list(range(1, 51))).all()
    merged.to_csv(OUT, index=False)
    old = base.set_index(["cancer", "rank"])["gene"]
    new = merged.set_index(["cancer", "rank"])["gene"]
    diff = [(c, r, old[(c, r)], new[(c, r)]) for (c, r) in old.index if old[(c, r)] != new[(c, r)]]
    print(f"Wrote {OUT.relative_to(ROOT)}; positions changed by the vocabulary correction: {len(diff)}")
    for c, r, a, b in diff:
        print(f"  {c} rank {r}: {a} -> {b}")


if __name__ == "__main__":
    main()
