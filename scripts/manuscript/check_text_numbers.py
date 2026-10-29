"""
check_text_numbers.py -- confirm that the headline numbers in the manuscript text equal the committed outputs.

Each check recomputes a value from a committed output file and requires its formatted string to appear in the
manuscript (abstract, results, supplement or tables). It does not prove every number in the paper is right; it
catches the failure that matters most, a number left over from a different build or an earlier output.
Exit code 1 if any check fails.

Run: python -X utf8 scripts/manuscript/check_text_numbers.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "analysis"))
import interaction_all19 as ia  # noqa: E402

TEXT = (ROOT / "manuscript" / "OncoDigger_manuscript1.1.md").read_text(encoding="utf-8")
INT = ROOT / "data" / "derived" / "intogen_rescoring"
fails: list[str] = []


def need(label: str, *needles: str) -> None:
    missing = [n for n in needles if n not in TEXT]
    print(("PASS  " if not missing else "FAIL  ") + label + (f"   missing: {missing}" if missing else ""))
    if missing:
        fails.append(label)


def main() -> None:
    G, E = ia.load_intogen()
    man = __import__("json").loads((ROOT / "data/provenance/corpus_build_2026-09-27.json").read_text(encoding="utf-8"))
    can = pd.read_csv(ROOT / "data/canonical/precision_at_k_detail.csv")
    t = can[can["rank"] <= 10]
    anyh = int(t.gene.isin(G).sum())
    exh = int(sum((g, c) in E for g, c in zip(t.gene, t.cancer)))
    need("corpus size", f"{man['total_documents']:,}")
    need("headline any-cancer", f"{anyh}/190", f"{100 * anyh / 190:.1f}%")
    need("headline exact", f"{exh}/190", f"{100 * exh / 190:.1f}%")

    st = pd.read_csv(INT / "statistical_support_reference_summary.csv")
    for met, lab in (("any_cancer_P@10", "statistical any"), ("exact_cancer_concordance", "statistical exact")):
        r = st[(st.system == "OncoDigger") & (st.reference == "statistical") & (st.metric == met)].iloc[0]
        need(lab, f"{int(r.hits)}/190 = {r.percent:.1f}%")

    ci = pd.read_csv(INT / "interaction_all19_bootstrap_ci.csv")
    for _, r in ci.iterrows():
        need(f"interaction {r.metric}", f"+{r.interaction_pp:.1f}", f"{r.ci_lo_pp:.1f}", f"{r.ci_hi_pp:.1f}")

    bc = pd.read_csv(INT / "precision_at_k_bootstrap_ci_intogen.csv").set_index("k").loc[10]
    need("bootstrap CI for Precision@10", f"{100 * bc['ci_low_2.5']:.1f}-{100 * bc['ci_high_97.5']:.1f}%")

    # one canonical cancer bootstrap everywhere: the same estimate must carry the same interval in every place it appears
    sm = pd.read_csv(INT / "statistical_support_reference_summary.csv")
    r = sm[(sm.system == "OncoDigger") & (sm.reference == "published") & (sm.metric == "exact_cancer_concordance")].iloc[0]
    need("exact concordance CI", f"95% CI {100 * r.boot_ci_lo:.1f}-{100 * r.boot_ci_hi:.1f}%")
    pp = pd.read_csv(INT / "statistical_support_reference_paired_vs_oncodigger.csv")
    for sysname in ("CancerMine", "ChatGPT (one response)", "Pooled cancer-blind control"):
        r = pp[(pp.system == sysname) & (pp.reference == "published") & (pp.metric == "exact_cancer_concordance")].iloc[0]
        need(f"paired exact CI vs {sysname}", f"{r.boot_ci_lo_pp:.1f}", f"{r.boot_ci_hi_pp:.1f}")
    r = pp[(pp.system == "CancerMine") & (pp.reference == "published") & (pp.metric == "exact_cancer_concordance")].iloc[0]
    need("CancerMine advantage CI, every place", f"(95% CI {r.boot_ci_lo_pp:.1f}-{r.boot_ci_hi_pp:.1f} pp,", f"CI. {r.boot_ci_lo_pp:.1f}–{r.boot_ci_hi_pp:.1f} pp")
    stale = [x for x in ("8.9-19.5", "8.9–19.5") if x in TEXT]
    print(("FAIL  " if stale else "PASS  ") + "no stale CancerMine interval 8.9-19.5")
    if stale:
        fails.append("stale CI")

    lw = pd.read_csv(INT / "intogen_literature_weight_audit.csv")
    ns = int(lw.engine_rank.notna().sum())
    need("surfaced / never surfaced", f"{ns}", f"{1740 - ns}")
    pres = pd.read_csv(INT / "intogen_full_corpus_presence_summary.csv").set_index("category")
    need("presence split", f"{int(pres.loc['present_but_outside_window'].n)}", f"{int(pres.loc['truly_absent_from_corpus'].n)}")

    mc = pd.read_csv(INT / "intogen_miss_root_cause_summary.csv")
    need("root-cause total", f"{int(mc.n.sum())} pairs")

    sec = __import__("json").loads((ROOT / "data/derived/secondary_stats.json").read_text(encoding="utf-8"))
    need("mean pairwise Jaccard", f"{sec['top10_pairwise_jaccard_mean']:.3f}")
    need("rank AUC", f"{sec['rank_auc']['mean']:.3f}", f"{sec['rank_auc']['min']:.3f}", f"{sec['rank_auc']['max']:.3f}")

    ms = pd.read_csv(INT / "vocab_correction_supporting_papers.csv")
    met = ms[(ms.gene == "MET")]
    need("MET supporting papers before/after", f"from {int(met.iloc[0].supporting_papers)} to {int(met.iloc[1].supporting_papers)}")
    strict = pd.read_csv(INT / "vocab_corrected_met_strict.csv").iloc[0]
    need("strict MET condition", f"{int(strict.met_postings_retained):,} postings", f"rank {int(strict.met_rank)}")

    pj = pd.read_csv(ROOT / "data/derived/pool_size_jaccard.csv").set_index("comparison")["mean_top10_jaccard"]
    need("pool-size Jaccard", f"{pj['250 vs 1000']:.3f}", f"{pj['750 vs 1000']:.3f}", f"{pj['2000 vs 1000']:.3f}", f"{pj['mean of adjacent']:.3f}")

    cm = pd.read_csv(INT / "cutoff_matched_summary.csv")
    for lab in ("OncoDigger (corpus to 2022)", "OncoDigger (corpus to 2023)"):
        r = cm[(cm.system == lab) & (cm.reference == "published") & (cm.metric == "exact_cancer_concordance")].iloc[0]
        need(f"cutoff exact {lab}", f"{int(r.hits)}/190 = {r.percent:.1f}%")

    tc = pd.read_csv(ROOT / "data/derived/vocab_rescoring/title_case_blind_spot_check.csv")
    for g in ("MET", "REST"):
        s = tc[tc.symbol == g]
        need(f"title-case check {g}", f"{int(s.retained_title_cased_title_only.sum())} of {int(s.retained.sum()):,}")

    # post hoc proximity-adjusted ranking (Supplementary Results S15): numbers come from data/derived/proximity/
    px = ROOT / "data" / "derived" / "proximity"
    sw = pd.read_csv(px / "summary_by_weight.csv", dtype=str).fillna("")
    f1 = lambda x: f"{float(x):.1f}"
    sgn = lambda x: f"{float(x):+.1f}".replace("-", "\u2212")
    w0, w3 = sw.iloc[0], sw[sw.weight == "3"].iloc[0]
    need("proximity w=3 exact", f"{w0.exact_10}/190 ({f1(w0.exact_pct)}%)", f"{w3.exact_10}/190 ({f1(w3.exact_pct)}%)",
         f"{sgn(w3.exact_diff_pp)} percentage points (95% cancer-bootstrap CI {sgn(w3.exact_ci_lo)} to {sgn(w3.exact_ci_hi)})")
    need("proximity w=3 any", f"{w0.any_10}/190 ({f1(w0.any_pct)}%)", f"{w3.any_10}/190 ({f1(w3.any_pct)}%)")
    for w in ("0.25", "1", "10"):
        rw = sw[sw.weight == w].iloc[0]
        rt6 = (ROOT / "data" / "supplementary_tables" / "RT6_formerS21_proximity_by_weight.csv").read_text(encoding="utf-8")
        ok_row = f"\n{w},{rw.exact_10} ({f1(rw.exact_pct)}%)" in rt6
        print(("PASS  " if ok_row else "FAIL  ") + f"proximity weight {w} row of Repository Table RT6")
        if not ok_row:
            fails.append(f"proximity weight {w} RT6 row")
    sel = pd.read_csv(px / "selection_honest.csv").iloc[0]
    need("proximity selection-aware", f"{sgn(sel.leave_one_cancer_out_gain_pp)} points",
         f"{sgn(sel.oob_bootstrap_gain_pp)} points (95% interval {sgn(sel.oob_ci_lo)} to {sgn(sel.oob_ci_hi)})")
    mc = pd.read_csv(px / "vs_matched_counting.csv").set_index("metric")
    ex, an = mc.loc["intogen_exact"], mc.loc["intogen_any"]
    asg = lambda x: f"{float(x):+.1f}"      # the manuscript prints ASCII hyphen-minus in this paragraph
    need("proximity vs matched counting", f"{asg(ex.diff_pp)} points on exact concordance (95% cancer-bootstrap CI {asg(ex.ci_lo)} to {asg(ex.ci_hi)}, "
         f"with {int(ex.cancers_better)} cancers better, {int(ex.cancers_equal)} equal and {int(ex.cancers_worse)} worse)",
         f"{asg(an.diff_pp)} points on any-cancer concordance (CI {asg(an.ci_lo)} to {asg(an.ci_hi)})")
    # control: the w = 0 ranking used by the proximity analysis equals the main (canonical) top 10 at every slot
    r0 = pd.read_csv(px / "ranking_w0.csv")
    m0 = can[can["rank"] <= 10].merge(r0[r0["rank"] <= 10], on=["cancer", "rank"], suffixes=("_c", "_p"))
    same = int((m0.gene_c.str.upper() == m0.gene_p.str.upper()).sum())
    print(("PASS  " if same == 190 else "FAIL  ") + f"proximity w=0 control equals canonical top 10 ({same}/190)")
    if same != 190:
        fails.append("proximity control")

    print(f"\n{len(fails)} check(s) failed" if fails else "\nAll text-number checks passed")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
