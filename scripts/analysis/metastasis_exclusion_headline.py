"""
metastasis_exclusion_headline.py -- the paper's headline numbers for all 19 cancers computed on the main ranking and on the ranking without
/secondary-indexed records (post hoc sensitivity analysis, Supplementary Results S16). Same IntOGen reference, mapping and bootstrap (19 cancers
resampled, 2,000 draws, seed 42) as the main analysis. Comparators (CancerMine, DISEASES, PubTator, ChatGPT) and the pooled control are external or
unchanged and were not recomputed.

Inputs : data/derived/query_planes/mutation.csv, data/derived/metastatic_exclusion/ranking*/mutation.csv, data/derived/pooled_control_mutation.csv
Output : data/derived/metastatic_exclusion/headline_main_vs_excluded.csv
Run    : python -X utf8 scripts/analysis/metastasis_exclusion_headline.py
"""
import importlib.util
import itertools
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "data" / "derived" / "metastatic_exclusion"
spec = importlib.util.spec_from_file_location("ir", ROOT / "scripts" / "analysis" / "intogen_full_rescoring.py")
ir = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ir)
G, E = ir.load_intogen()
G = {g.upper() for g in G}
E = {(g.upper(), c) for g, c in E}
CANCERS = ir.STUDY_ORDER


def rank_lists(df):
    return {c: g.sort_values("rank").gene.tolist() for c, g in df.groupby("cancer")}


def per_cancer(lists, k, kind):
    f = (lambda g, c: g.upper() in G) if kind == "any" else (lambda g, c: (g.upper(), c) in E)
    return np.array([np.mean([f(g, c) for g in lists[c][:k]]) for c in CANCERS])


def headline(lists, pooled10):
    out = {}
    any10, ex10 = per_cancer(lists, 10, "any"), per_cancer(lists, 10, "exact")
    for name, arr in (("any", any10), ("exact", ex10)):
        m, lo, hi = ir.case_bootstrap_ci(arr)
        out[f"{name}_hits_of_190"] = int(round(arr.sum() * 10))
        out[f"{name}_pct"] = round(100 * m, 1)
        out[f"{name}_ci"] = f"{100 * lo:.1f}-{100 * hi:.1f}"
    out["cancers_with_10_of_10_any"] = int((any10 == 1).sum())
    out["lowest_exact_cancers"] = ", ".join(c for c, v in zip(CANCERS, ex10) if v == ex10.min()) + f" ({int(round(ex10.min() * 10))}/10)"
    for k in (5, 10, 15, 20, 30, 40, 50):
        out[f"mean_P@{k}"] = round(100 * per_cancer(lists, k, "any").mean(), 1)
    pairs_total = len(E)
    for k in (10, 20, 50):
        out[f"pooled_recall@{k}"] = f"{sum(sum((g.upper(), c) in E for g in lists[c][:k]) for c in CANCERS)}/{pairs_total}"
    t = {c: set(lists[c][:10]) for c in CANCERS}
    out["mean_pairwise_cancer_top10_jaccard"] = round(float(np.mean([len(t[a] & t[b]) / len(t[a] | t[b]) for a, b in itertools.combinations(CANCERS, 2)])), 3)
    extra = {c: [g for g in lists[c][:10] if g not in pooled10] for c in CANCERS}
    out["extra_genes_vs_pooled_per_cancer"] = round(float(np.mean([len(v) for v in extra.values()])), 2)
    out["extra_genes_IntOGen_any_per_cancer"] = round(float(np.mean([sum(g.upper() in G for g in v) for v in extra.values()])), 2)
    n_extra = sum(len(v) for v in extra.values())
    n_extra_exact = sum((g.upper(), c) in E for c, v in extra.items() for g in v)
    out["extra_genes_exact"] = f"{n_extra_exact}/{n_extra} ({100 * n_extra_exact / n_extra:.1f}%)"
    nonmatch = [(c, g) for c in CANCERS for g in lists[c][:10] if (g.upper(), c) not in E]
    other = sum(g.upper() in G for _, g in nonmatch)
    out["nonmatching_top10_pairs"] = f"{len(nonmatch)} (driver in another cancer {other}, not a driver anywhere {len(nonmatch) - other})"
    return out, any10, ex10


def main() -> None:
    base = pd.read_csv(ROOT / "data" / "derived" / "query_planes" / "mutation.csv")
    new = pd.concat([pd.read_csv(D / "ranking" / "mutation.csv"), pd.read_csv(D / "ranking_other14" / "mutation.csv")], ignore_index=True)
    assert new.cancer.nunique() == 19 and len(new) == 950
    pooled10 = set(pd.read_csv(ROOT / "data" / "derived" / "pooled_control_mutation.csv").sort_values("rank").gene.tolist()[:10])
    hb, ab, eb = headline(rank_lists(base), pooled10)
    hn, an, en = headline(rank_lists(new), pooled10)
    rows = [{"measure": k, "main_ranking": hb[k], "without_secondary_records": hn[k]} for k in hb]
    for name, a, b in (("any", ab, an), ("exact", eb, en)):
        m, lo, hi, _ = ir.paired_bootstrap_ci(b, a)
        rows.append({"measure": f"paired change {name}, excluded minus main (points, 95% CI; cancers better/equal/worse)",
                     "main_ranking": "", "without_secondary_records": f"{100 * m:+.1f} ({100 * lo:+.1f} to {100 * hi:+.1f}); {int((b > a).sum())}/{int((b == a).sum())}/{int((b < a).sum())}"})
    df = pd.DataFrame(rows)
    df.to_csv(D / "headline_main_vs_excluded.csv", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
