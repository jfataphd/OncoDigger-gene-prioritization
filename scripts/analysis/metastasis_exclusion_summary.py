"""
metastasis_exclusion_summary.py -- compare the main ranking with the ranking computed without /secondary-indexed records, for all 19 cancers
(post hoc sensitivity analysis, Supplementary Results S16). Nothing in the main analysis is changed.

The sensitivity ranking is stored in two files because it was run in two batches with the same pipeline and settings:
  data/derived/metastatic_exclusion/ranking/mutation.csv          Liver, Brain, Lung, Kidney, Melanoma
  data/derived/metastatic_exclusion/ranking_other14/mutation.csv  the other 14 cancers
Main ranking: data/derived/query_planes/mutation.csv (frozen build, policy v5; its top-50 lists equal those of mutation_vocab_corrected.csv).

Outputs (data/derived/metastatic_exclusion/): comparison_per_cancer.csv, gene_movements_top10.csv, liver_selected_genes.csv, comparison_overall.csv
Run    : python -X utf8 scripts/analysis/metastasis_exclusion_summary.py
"""
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "analysis"))
import interaction_all19 as ia  # noqa: E402

D = ROOT / "data" / "derived" / "metastatic_exclusion"
LABELS = {"breast": "Breast Cancer", "lung": "Lung Cancer", "colorectal": "Colorectal Cancer", "prostate": "Prostate Cancer", "melanoma": "Melanoma",
          "bladder": "Bladder Cancer", "kidney": "Kidney Cancer", "pancreas": "Pancreatic Cancer", "liver": "Liver Cancer",
          "stomach": "Stomach Cancer", "esophageal": "Esophageal Cancer", "ovarian": "Ovarian Cancer", "endometrial": "Endometrial Cancer",
          "cervical": "Cervical Cancer", "thyroid": "Thyroid Cancer", "brain": "Brain Cancer", "leukemia": "Leukemia", "lymphoma": "Lymphoma",
          "myeloma": "Myeloma"}
G, E = ia.load_intogen()


def lists(df, cancer, k):
    return df[df.cancer == cancer].sort_values("rank").gene.tolist()[:k]


def jacc(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b)


def main() -> None:
    base = pd.read_csv(ROOT / "data" / "derived" / "query_planes" / "mutation.csv")
    new = pd.concat([pd.read_csv(D / "ranking" / "mutation.csv"), pd.read_csv(D / "ranking_other14" / "mutation.csv")], ignore_index=True)
    assert new.cancer.nunique() == 19 and len(new) == 950, (new.cancer.nunique(), len(new))
    summ = pd.read_csv(D / "secondary_summary.csv").set_index("collection")
    rows, moves = [], []
    for key, label in LABELS.items():
        b10, n10, b50, n50 = lists(base, label, 10), lists(new, label, 10), lists(base, label, 50), lists(new, label, 50)
        hit = lambda l: sum(g.upper() in G for g in l)
        ex = lambda l: sum((g.upper(), label) in E for g in l)
        rows.append({
            "cancer": label, "removed_records": int(summ.loc[key, "in_frozen_collection"]),
            "removed_pct_of_frozen_collection": float(summ.loc[key, "share_of_frozen_collection_pct"]),
            "top10_genes_retained": len(set(b10) & set(n10)), "top10_jaccard": round(jacc(b10, n10), 3), "top50_jaccard": round(jacc(b50, n50), 3),
            "top10_same_order": b10 == n10, "top50_identical": b50 == n50,
            "any_main": hit(b10), "any_excluded": hit(n10), "exact_main": ex(b10), "exact_excluded": ex(n10),
        })
        for g in b10:
            if g not in n10:
                moves.append({"cancer": label, "gene": g, "move": "leaves top 10", "main_rank": b50.index(g) + 1 if g in b50 else None,
                              "excluded_rank": n50.index(g) + 1 if g in n50 else None})
        for g in n10:
            if g not in b10:
                moves.append({"cancer": label, "gene": g, "move": "enters top 10", "main_rank": b50.index(g) + 1 if g in b50 else None,
                              "excluded_rank": n50.index(g) + 1})
    r = pd.DataFrame(rows)
    r.to_csv(D / "comparison_per_cancer.csv", index=False)
    pd.DataFrame(moves).to_csv(D / "gene_movements_top10.csv", index=False)
    sel = []
    b50, n50 = lists(base, "Liver Cancer", 50), lists(new, "Liver Cancer", 50)
    for g in ("KRAS", "BRAF", "NRAS", "APC", "TP53", "CTNNB1", "TERT", "EGFR", "PIK3CA", "HFE"):
        sel.append({"gene": g, "main_rank": b50.index(g) + 1 if g in b50 else None, "excluded_rank": n50.index(g) + 1 if g in n50 else None})
    pd.DataFrame(sel).to_csv(D / "liver_selected_genes.csv", index=False)

    def mean_pair_jacc(df):
        t = {c: set(lists(df, c, 10)) for c in LABELS.values()}
        return float(np.mean([len(t[a] & t[b]) / len(t[a] | t[b]) for a, b in itertools.combinations(t, 2)]))

    def pk(df, k):
        return float(np.mean([np.mean([g.upper() in G for g in lists(df, c, k)]) for c in LABELS.values()]) * 100)

    overall = {
        "cancers_with_changed_top10_set": int((r.top10_genes_retained < 10).sum()),
        "cancers_with_identical_top50": int(r.top50_identical.sum()),
        "top10_slots_changed": int((10 - r.top10_genes_retained).sum()), "of_slots": 190,
        "any_main": int(r.any_main.sum()), "any_excluded": int(r.any_excluded.sum()),
        "exact_main": int(r.exact_main.sum()), "exact_excluded": int(r.exact_excluded.sum()),
        "mean_top10_jaccard_main_vs_excluded": round(float(r.top10_jaccard.mean()), 3),
        "mean_pairwise_cancer_top10_jaccard_main": round(mean_pair_jacc(base), 3),
        "mean_pairwise_cancer_top10_jaccard_excluded": round(mean_pair_jacc(new), 3),
        "mean_P@20_main": round(pk(base, 20), 1), "mean_P@20_excluded": round(pk(new, 20), 1),
        "mean_P@50_main": round(pk(base, 50), 1), "mean_P@50_excluded": round(pk(new, 50), 1),
    }
    pd.DataFrame([overall]).T.rename(columns={0: "value"}).to_csv(D / "comparison_overall.csv")
    print(r.to_string(index=False))
    print("\nOVERALL:", overall)
    print("\ntop-10 movements:\n", pd.DataFrame(moves).to_string(index=False))


if __name__ == "__main__":
    main()
