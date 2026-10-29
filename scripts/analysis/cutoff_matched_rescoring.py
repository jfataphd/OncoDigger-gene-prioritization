"""
cutoff_matched_rescoring.py -- literature-cutoff-matched sensitivity analysis.

Question: how much of OncoDigger's advantage over CancerMine reflects access to newer literature?
CancerMine's frozen snapshot (Zenodo 7689627, 2023-03-01) cannot contain later papers, while
OncoDigger's corpora run to September 2026. Here OncoDigger is re-run, with the algorithm,
vocabulary and settings unchanged, on corpora truncated by publication year, and the rankings are
scored against the SAME fixed IntOGen reference (release 2024.09.20) as everything else.

The rankings in data/derived/cutoff_matched_rankings/ (the *_vocab_corrected.csv files: the same MET/REST/FH
redaction as the main analysis, applied by cutoff_matched_vocab_correction.py to the truncated Lung, Thyroid and
Kidney collections; the uncorrected files are kept and give identical counts) were produced by scripts/analysis/
run_query_planes.py (unchanged) on corpora built by scripts/analysis/build_date_truncated_corpora.py,
which recomputes the corpus size, IDF, average document length and enrichment background counts on
the truncated set (so this is a true "corpus as of year X" run, not a retrieval-only filter).
  cutoff 2022: publication_year <= 2022 (slightly stricter than CancerMine's March 2023 snapshot)
  cutoff 2023: publication_year <= 2023 (slightly looser, includes ~10 months CancerMine lacks)
The two cutoffs bracket the snapshot date. The IntOGen reference is not truncated, so both
OncoDigger and CancerMine are scored against a reference that contains post-cutoff evidence.

Reference definitions and the cancer-bootstrap scheme (2000 resamples of the 19 cancers, seed 42)
are imported unchanged from intogen_statistical_support_rescoring.py.

Run: python -X utf8 scripts/analysis/cutoff_matched_rescoring.py
Outputs (data/derived/intogen_rescoring/):
  cutoff_matched_summary.csv, cutoff_matched_paired.csv, cutoff_matched_overlap.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import intogen_statistical_support_rescoring as base  # noqa: E402

ROOT = base.ROOT
RANK_DIR = ROOT / "data" / "derived" / "cutoff_matched_rankings"
CUTS = {"OncoDigger (corpus to 2022)": RANK_DIR / "mutation_corpus_to_2022_vocab_corrected.csv",
        "OncoDigger (corpus to 2023)": RANK_DIR / "mutation_corpus_to_2023_vocab_corrected.csv"}
PAPER = "OncoDigger (paper, corpus to 2026)"
COMPARATORS = [s for s in base.SYSTEMS if s != "OncoDigger"]
METRICS = (("any", "any_cancer_P@10"), ("exact", "exact_cancer_concordance"))


def load_cut(path: Path) -> pd.DataFrame:
    t = pd.read_csv(path)
    t = t[t["rank"] <= 10][["cancer", "rank", "gene"]].copy()
    t["gene"] = t["gene"].astype(str).str.strip().str.upper()
    assert len(t) == 190 and t.cancer.nunique() == 19, path
    return t


def per_cancer(t: pd.DataFrame, any_genes, exact) -> dict:
    t = t.copy()
    t["any"] = t["gene"].isin(any_genes)
    t["exact"] = [(g, c) in exact for g, c in zip(t["gene"], t["cancer"])]
    out = {}
    for m, _ in METRICS:
        out[m] = (t.groupby("cancer")[m].mean().reindex(base.CANCERS).to_numpy(), int(t[m].sum()))
    return out


def main():
    refs = base.load_references()
    systems = {PAPER: base.load_system("canon", None)}
    systems.update({k: load_cut(v) for k, v in CUTS.items()})
    for name in COMPARATORS:
        kind, fname = base.SYSTEMS[name]
        systems[name] = base.load_system(kind, fname)

    scored = {}
    summary = []
    for sname, t in systems.items():
        for rname, (anyg, exact, _) in refs.items():
            pc = per_cancer(t, anyg, exact)
            for m, label in METRICS:
                v, hits = pc[m]
                lo, hi = base.boot_ci(v)
                scored[(sname, rname, m)] = v
                summary.append(dict(system=sname, reference=rname, metric=label, hits=hits, n=190,
                                    percent=round(100 * hits / 190, 1), boot_ci_lo=lo, boot_ci_hi=hi))
    summary = pd.DataFrame(summary)

    paired = []
    for cut in list(CUTS) + [PAPER]:
        for comp in COMPARATORS:
            for rname in refs:
                for m, label in METRICS:
                    d = scored[(cut, rname, m)] - scored[(comp, rname, m)]
                    lo, hi = base.boot_ci(d)
                    paired.append(dict(oncodigger=cut, comparator=comp, reference=rname, metric=label,
                                       mean_diff_pp=round(100 * float(d.mean()), 1),
                                       boot_ci_lo_pp=round(100 * lo, 1), boot_ci_hi_pp=round(100 * hi, 1),
                                       cancers_oncodigger_better=int((d > 0).sum()), ties=int((d == 0).sum()),
                                       cancers_oncodigger_worse=int((d < 0).sum())))
    paired = pd.DataFrame(paired)

    paper = systems[PAPER]
    overlap = []
    for cut in CUTS:
        t = systems[cut]
        same_set = sum(len(set(paper[paper.cancer == c].gene) & set(t[t.cancer == c].gene)) for c in base.CANCERS)
        same_order = sum(paper[paper.cancer == c].sort_values("rank").gene.tolist()
                         == t[t.cancer == c].sort_values("rank").gene.tolist() for c in base.CANCERS)
        overlap.append(dict(oncodigger=cut, top10_slots_shared_with_paper=same_set, of=190,
                            cancers_with_identical_ordered_list=same_order, of_cancers=19))
    overlap = pd.DataFrame(overlap)

    out = ROOT / "data" / "derived" / "intogen_rescoring"
    summary.to_csv(out / "cutoff_matched_summary.csv", index=False)
    paired.to_csv(out / "cutoff_matched_paired.csv", index=False)
    overlap.to_csv(out / "cutoff_matched_overlap.csv", index=False)

    pd.set_option("display.width", 220)
    keep = [PAPER] + list(CUTS) + ["CancerMine", "ChatGPT (one response)", "Pooled cancer-blind control",
                                   "Matched frequency-only (clean vocabulary)"]
    for label in ("any_cancer_P@10", "exact_cancer_concordance"):
        s = summary[(summary.metric == label) & summary.system.isin(keep)]
        print(f"\n{label} (hits of 190)")
        print(s.pivot(index="system", columns="reference", values="hits").reindex(keep).to_string())
    print()
    print(overlap.to_string(index=False))
    show = paired[(paired.comparator == "CancerMine")]
    print("\nPaired vs CancerMine")
    print(show.to_string(index=False))


if __name__ == "__main__":
    main()
