"""
proximity_summary.py -- every number reported for the post hoc proximity-adjusted ranking (Supplementary Results S15).

Inputs (data/derived/proximity/, produced by proximity_boost_run.py on the September-27 build with the MET/REST/FH-corrected collections):
  ranking_w0.csv  ranking_w0.25.csv  ranking_w1.csv  ranking_w3.csv  ranking_w10.csv   cancer, rank, gene (top 50)
  ranking_rule_w0.csv  ranking_rule_w3.csv                                               alias-rule sensitivity
  scores_w3.csv                                                                          cancer, gene, supporting_papers, original_relevance, proximity, adjusted_relevance
Outputs (same folder): summary_by_weight.csv, summary_per_cancer.csv, gene_moves_w3.csv, selection_honest.csv, alias_rule_summary.csv,
proximity_by_status.csv.
All concordance figures use the same IntOGen reference, cancer mapping and 19-cancer case bootstrap (2,000 draws, seed 42) as the main analysis
(intogen_full_rescoring.py, intogen_statistical_support_rescoring.py). Run: python -X utf8 scripts/analysis/proximity_summary.py
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "data" / "derived" / "proximity"


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / "analysis" / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


ir = load("intogen_full_rescoring")
ss = load("intogen_statistical_support_rescoring")
GA, PAIRS = ir.load_intogen()
REFS = ss.load_references()
ORDER = ir.STUDY_ORDER
WEIGHTS = ["0.25", "1", "3", "10"]


def rk(path):
    d = pd.read_csv(path)
    d["gene"] = d.gene.astype(str).str.upper()
    return d


def top10(d):
    return d[d["rank"] <= 10][["cancer", "rank", "gene"]]


def per_cancer(d):
    t = ir.score(top10(d)[["cancer", "gene"]], GA, PAIRS)
    return t.groupby("cancer")[["intogen_any", "intogen_exact"]].sum().reindex(ORDER)


def stat_counts(d):
    t = top10(d)
    a, e, _ = REFS["statistical"]
    return int(t.gene.isin(a).sum()), int(sum((g, c) in e for g, c in zip(t.gene, t.cancer)))


def shared_top10(a, b):
    return int(sum(len(set(top10(a)[top10(a).cancer == c].gene) & set(top10(b)[top10(b).cancer == c].gene)) for c in ORDER))


def boot_diff(x, y):
    m, lo, hi, _ = ir.paired_bootstrap_ci(x.values / 10.0, y.values / 10.0)
    return 100 * m, 100 * lo, 100 * hi


def main():
    base = rk(D / "ranking_w0.csv")
    runs = {w: rk(D / f"ranking_w{w}.csv") for w in WEIGHTS}
    pcb = per_cancer(base)

    rows = []
    for name, d in [("0 (frozen ranking)", base)] + [(w, runs[w]) for w in WEIGHTS]:
        pc = per_cancer(d)
        sa, se = stat_counts(d)
        r = {"weight": name, "any_10": int(pc.intogen_any.sum()), "exact_10": int(pc.intogen_exact.sum()),
             "any_pct": round(100 * pc.intogen_any.sum() / 190, 1), "exact_pct": round(100 * pc.intogen_exact.sum() / 190, 1),
             "stat_any_10": sa, "stat_exact_10": se, "top10_genes_kept": shared_top10(base, d)}
        if name == "0 (frozen ranking)":
            r.update({"exact_diff_pp": "", "exact_ci_lo": "", "exact_ci_hi": "", "any_diff_pp": "", "any_ci_lo": "", "any_ci_hi": "",
                      "cancers_better_exact": "", "cancers_worse_exact": ""})
        else:
            m, lo, hi = boot_diff(pc.intogen_exact, pcb.intogen_exact)
            ma, loa, hia = boot_diff(pc.intogen_any, pcb.intogen_any)
            r.update({"exact_diff_pp": round(m, 1), "exact_ci_lo": round(lo, 1), "exact_ci_hi": round(hi, 1),
                      "any_diff_pp": round(ma, 1), "any_ci_lo": round(loa, 1), "any_ci_hi": round(hia, 1),
                      "cancers_better_exact": int((pc.intogen_exact > pcb.intogen_exact).sum()),
                      "cancers_worse_exact": int((pc.intogen_exact < pcb.intogen_exact).sum())})
        rows.append(r)
    pd.DataFrame(rows).to_csv(D / "summary_by_weight.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))

    # per-cancer table at w = 3
    p3 = per_cancer(runs["3"])
    pcdf = pd.DataFrame({"cancer": ORDER, "exact_frozen": pcb.intogen_exact.values.astype(int), "exact_w3": p3.intogen_exact.values.astype(int),
                         "any_frozen": pcb.intogen_any.values.astype(int), "any_w3": p3.intogen_any.values.astype(int)})
    pcdf["exact_change"] = pcdf.exact_w3 - pcdf.exact_frozen
    pcdf["top10_genes_kept"] = [len(set(top10(base)[top10(base).cancer == c].gene) & set(top10(runs["3"])[top10(runs["3"]).cancer == c].gene)) for c in ORDER]
    pcdf.to_csv(D / "summary_per_cancer.csv", index=False)

    # genes that enter / leave the top 10 at w = 3
    moves = []
    for c in ORDER:
        b = set(top10(base)[top10(base).cancer == c].gene)
        n = set(top10(runs["3"])[top10(runs["3"]).cancer == c].gene)
        for g in sorted(b - n):
            moves.append((c, g, "leaves", g in GA, (g, c) in PAIRS))
        for g in sorted(n - b):
            moves.append((c, g, "enters", g in GA, (g, c) in PAIRS))
    mv = pd.DataFrame(moves, columns=["cancer", "gene", "move", "intogen_any", "intogen_exact"])
    mv.to_csv(D / "gene_moves_w3.csv", index=False)
    print("\ngene moves at w=3:", len(mv), "(", int((mv.move == 'leaves').sum()), "leave )")
    print(mv.groupby("move")[["intogen_any", "intogen_exact"]].sum().to_string())

    # selection-honest summaries over the four weights (exact concordance)
    names = WEIGHTS
    M = np.array([(per_cancer(runs[w]).intogen_exact.values - pcb.intogen_exact.values) / 10.0 for w in names])    # settings x cancers, in proportion units
    loco_total, picks = 0.0, []
    for i in range(19):
        rest = [j for j in range(19) if j != i]
        k = int(np.argmax(M[:, rest].sum(axis=1)))
        picks.append(names[k])
        loco_total += M[k, i]
    rng = np.random.default_rng(2026)
    oob, inbag = [], []
    for _ in range(5000):
        idx = rng.integers(0, 19, 19)
        inb = np.zeros(19, bool)
        inb[idx] = True
        if inb.all():
            continue
        sc = M[:, idx].mean(axis=1)
        k = int(np.argmax(sc))
        oob.append(M[k][~inb].mean())
        inbag.append(sc.max())
    oob = np.array(oob)
    sel = pd.DataFrame([{
        "leave_one_cancer_out_gain_pp": round(100 * loco_total / 19, 1), "loco_weights_chosen": str(pd.Series(picks).value_counts().to_dict()),
        "oob_bootstrap_gain_pp": round(100 * oob.mean(), 1), "oob_ci_lo": round(100 * np.percentile(oob, 2.5), 1),
        "oob_ci_hi": round(100 * np.percentile(oob, 97.5), 1), "oob_share_positive_pct": round(100 * (oob > 0).mean()),
        "optimism_pp": round(100 * (np.mean(inbag) - oob.mean()), 1)}])
    sel.to_csv(D / "selection_honest.csv", index=False)
    print("\n", sel.to_string(index=False))

    # alias-rule sensitivity
    rb, rs = rk(D / "ranking_rule_w0.csv"), rk(D / "ranking_rule_w3.csv")
    al = []
    for name, d, ref in (("frozen vocabulary, w = 0", base, None), ("frozen vocabulary, w = 3", runs["3"], base),
                         ("alias rule, w = 0", rb, base), ("alias rule, w = 3", rs, rb)):
        pc = per_cancer(d)
        r = {"setting": name, "any_10": int(pc.intogen_any.sum()), "exact_10": int(pc.intogen_exact.sum()),
             "stat_any_10": stat_counts(d)[0], "stat_exact_10": stat_counts(d)[1]}
        if ref is not None:
            pr = per_cancer(ref)
            m, lo, hi = boot_diff(pc.intogen_exact, pr.intogen_exact)
            r.update({"exact_diff_vs_reference_pp": round(m, 1), "ci_lo": round(lo, 1), "ci_hi": round(hi, 1)})
        al.append(r)
    mm = base.merge(rb, on=["cancer", "rank"])
    extra = {"setting": "alias rule vs frozen: slots changed (top 10 / top 50)",
             "any_10": int((mm[mm["rank"] <= 10].gene_x != mm[mm["rank"] <= 10].gene_y).sum()), "exact_10": int((mm.gene_x != mm.gene_y).sum())}
    al.append(extra)
    pd.DataFrame(al).to_csv(D / "alias_rule_summary.csv", index=False)
    print("\n", pd.DataFrame(al).to_string(index=False))

    # proximity score by IntOGen status (top-10 genes of the frozen ranking) and by rank band
    sc = pd.read_csv(D / "scores_w3.csv")
    sc["gene"] = sc.gene.astype(str).str.upper()
    t = top10(base).merge(sc[["cancer", "gene", "supporting_papers", "proximity"]], on=["cancer", "gene"], how="left")
    t["status"] = ["driver in this cancer" if (g, c) in PAIRS else "driver in another cancer" if g in GA else "not an IntOGen driver"
                   for g, c in zip(t.gene, t.cancer)]
    ps = t.groupby("status").agg(pairs=("gene", "size"), median_proximity=("proximity", "median"), median_supporting_papers=("supporting_papers", "median")).round(3)
    ps.to_csv(D / "proximity_by_status.csv")
    allr = base[base["rank"] <= 50].merge(sc[["cancer", "gene", "proximity"]], on=["cancer", "gene"], how="left")
    allr["band"] = np.where(allr["rank"] <= 10, "ranks 1-10", np.where(allr["rank"] <= 25, "ranks 11-25", "ranks 26-50"))
    print("\n", ps.to_string())
    print(allr.groupby("band").proximity.median().round(3).to_string())


if __name__ == "__main__":
    main()
