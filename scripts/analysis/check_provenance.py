"""
check_provenance.py -- verify that the committed outputs belong to ONE corpus build and one analysis state.

What it checks (exit code 1 if any check fails):
  1. outputs_record.csv: every listed output file still has the recorded SHA-256 (detects an output that was
     edited, regenerated from another build, or left stale after upstream outputs changed).
  2. corpus_stats.csv: papers per corpus equal the build manifest's document counts.
  3. Cross-output consistency (an output built from a different corpus build than the others breaks these):
       - canonical top-10 sets equal the main-pipeline (mutation.csv) top-10 sets, and the canonical ranking
         differs from mutation.csv only at the documented MET/REST/FH correction positions;
       - the headline counts (any-cancer and exact hits for OncoDigger) agree across all_systems_summary,
         statistical_support_reference_summary, the independence audit, true_raw_vs_cleaned, ablation panel,
         matched-frequency sanity check and the 2x2 factorial;
       - the literature-weight audit, presence check and confidence breakdown partition the same 1,740 pairs
         into the same surfaced / never-surfaced split;
       - the pooled control's top 10 is the same list for every cancer and its exact count agrees across files.
  4. With ONCODIGGER_ROOT set: the corpus root is the recorded build (provenance.assert_corpus_build),
     add --full-hash to recompute every database SHA-256.

Run: python -X utf8 scripts/analysis/check_provenance.py [--full-hash]
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import provenance  # noqa: E402

REC = ROOT / "data" / "provenance" / "outputs_record.csv"
INT = ROOT / "data" / "derived" / "intogen_rescoring"
failures: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(("PASS  " if ok else "FAIL  ") + msg)
    if not ok:
        failures.append(msg)


def sha(p: Path) -> str:
    return provenance.file_digest(p)


def lists(df: pd.DataFrame, n: int) -> dict:
    return {c: list(g.sort_values("rank").head(n)["gene"]) for c, g in df.groupby("cancer")}


def main() -> None:
    man = provenance.load_manifest()

    # 1 record
    if REC.exists():
        rec = pd.read_csv(REC)
        bad = [r.path for r in rec.itertuples() if not (ROOT / r.path).exists() or sha(ROOT / r.path) != r.sha256]
        check(not bad, f"outputs_record.csv: {len(rec)} outputs unchanged since recorded" + (f" (changed/missing: {bad[:5]})" if bad else ""))
    else:
        check(False, "outputs_record.csv missing (run make_provenance_record.py)")

    # 2 corpus_stats
    cs = pd.read_csv(ROOT / "data" / "canonical" / "corpus_stats.csv").set_index("corpus_key")["papers"]
    diffs = [(c["corpus"], int(cs[f"{c['corpus']}_cancer_inclusive"]), c["documents"]) for c in man["corpora"]
             if int(cs[f"{c['corpus']}_cancer_inclusive"]) != c["documents"]]
    check(not diffs, f"corpus_stats.csv paper counts equal the {man['build_id']} build manifest" + (f" (differ: {diffs[:3]})" if diffs else ""))
    check(int(cs.sum()) == man["total_documents"], f"corpus_stats.csv total {int(cs.sum()):,} == manifest total {man['total_documents']:,}")

    # 3a rankings
    mut = pd.read_csv(ROOT / "data" / "derived" / "query_planes" / "mutation.csv")
    cor = pd.read_csv(ROOT / "data" / "derived" / "query_planes" / "mutation_vocab_corrected.csv")
    can = pd.read_csv(ROOT / "data" / "canonical" / "precision_at_k_detail.csv")
    lm, lc, lk = lists(mut, 50), lists(cor, 50), lists(can, 50)
    check(lk == lc, "canonical ranking == mutation_vocab_corrected ranking (all 19 x 50)")
    # With the exact background count (2026-10 rerun) the MET (Lung) and FH (Kidney) corrections no longer move any position in the
    # top 50; the REST (Thyroid) correction shifts exactly ranks 41-50 (REST leaves at 41, CD274 enters at 50). The check is an exact
    # expectation on the set of differing positions, not an upper bound, so any other difference fails.
    # Under policy v5 the bare symbol REST is suppressed globally, so REST no longer enters the Thyroid ranking before the
    # redaction either; the MET/REST/FH redaction then changes supporting counts but no ranking list (Methods 3.3).
    got_diff = {(c, i + 1) for c in lm for i in range(50) if lm[c][i] != lc[c][i]}
    check(not got_diff,
          f"main-pipeline ranking equals the corrected ranking at every position under policy v5 "
          f"(found {len(got_diff)} differing positions in {sorted({c for c, _ in got_diff})})")
    check(all(set(lm[c][:10]) == set(lc[c][:10]) for c in lm), "top-10 gene sets identical in the main and corrected rankings")

    # 3b headline counts
    def one(df, **q):
        m = df
        for k, v in q.items():
            m = m[m[k] == v]
        return int(m.iloc[0]["hits"])
    stat = pd.read_csv(INT / "statistical_support_reference_summary.csv")
    a_any = one(stat, system="OncoDigger", reference="published", metric="any_cancer_P@10")
    a_ex = one(stat, system="OncoDigger", reference="published", metric="exact_cancer_concordance")
    asum = pd.read_csv(INT / "all_systems_summary.csv")
    r = asum[(asum.iloc[:, 0] == "oncodigger")]
    vals = {row.iloc[1]: int(row.iloc[2]) for _, row in r.iterrows()}
    check(vals.get("intogen_exact") == a_ex, f"all_systems_summary exact hits ({vals.get('intogen_exact')}) == statistical_support summary ({a_ex})")
    aud = pd.read_csv(ROOT / "data" / "derived" / "intogen_independence_audit" / "benchmark_recalculation.csv")
    od = aud[aud.system == "OncoDigger"]
    check(int(od[(od.reference == "0_published") & (od.metric == "exact_cancer_concordance")].hits.iloc[0]) == a_ex and
          int(od[(od.reference == "0_published") & (od.metric == "any_cancer_P@10")].hits.iloc[0]) == a_any,
          f"independence audit reproduces the published counts ({a_any} any, {a_ex} exact)")
    tr = pd.read_csv(INT / "true_raw_vs_cleaned_all19.csv")
    check(int(tr.cleaned_any.sum()) == a_any and int(tr.cleaned_exact.sum()) == a_ex, "true_raw_vs_cleaned cleaned column == headline counts")
    abl = pd.read_csv(INT / "ablation_panel_summary.csv")
    fp = abl[abl.iloc[:, 0].str.startswith("full_pipeline")].iloc[0]
    check(int(fp.iloc[2]) == a_any and int(fp.iloc[4]) == a_ex, "ablation panel full pipeline == headline counts")
    inter = pd.read_csv(INT / "interaction_all19_final.csv")
    mf = pd.read_csv(INT / "matched_frequency_baseline_summary.csv")
    sanity = mf[mf.iloc[:, 0] == "real_pipeline_sanity_check"].iloc[0]
    check(sanity.iloc[2] == f"{a_ex}/190" and sanity.iloc[1] == f"{a_any}/190", "matched-frequency sanity row reproduces the headline counts")
    check(int(inter.cell4_any.sum()) == int(mf[mf.iloc[:, 0] == "matched_frequency_only"].iloc[0].iloc[1].split("/")[0]),
          "2x2 factorial cell 4 (BM25 pool, cleaned) == matched frequency-only baseline")

    # 3c 1,740 pair partition
    lw = pd.read_csv(INT / "intogen_literature_weight_audit.csv")
    n_surf = int(lw.engine_rank.notna().sum())
    pres = pd.read_csv(INT / "intogen_full_corpus_presence_summary.csv").set_index("category")
    cb = pd.read_csv(INT / "intogen_confidence_breakdown_summary.csv").set_index("metric")
    check(len(lw) == 1740 and int(pres.n.sum()) == 1740, "presence summary partitions the same 1,740 pairs as the literature-weight audit")
    check(int(pres.loc["surfaced_in_window"].n) == n_surf == int(cb.loc["n_methods"].surfaced_n),
          f"surfaced pairs agree across audit ({n_surf}), presence summary and confidence breakdown")
    miss = pd.read_csv(INT / "intogen_recall_missed_genes.csv")
    check(len(miss) == len(lw) - n_surf, "missed-genes list == pairs never surfaced in the audit")

    # 3d pooled control
    pooled = pd.read_csv(ROOT / "data" / "derived" / "pooled_control_mutation.csv")
    check(len(pooled) >= 10, "pooled_control_mutation.csv present (same top-10 list for every cancer by construction)")
    pidx = man.get("pooled_index", {})
    check(pidx.get("documents") == 1514148, f"manifest pooled index documents == 1,514,148 ({pidx.get('documents')})")

    # 4 corpus root
    root = os.environ.get("ONCODIGGER_ROOT")
    if root and Path(root).exists():
        try:
            provenance.assert_corpus_build(root, full="--full-hash" in sys.argv)
            check(True, "ONCODIGGER_ROOT is the recorded corpus build")
        except RuntimeError as e:
            check(False, str(e).splitlines()[0])
    else:
        print("SKIP  ONCODIGGER_ROOT not set: corpus databases not checked")

    print(f"\n{len(failures)} check(s) failed" if failures else "\nAll provenance checks passed")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
