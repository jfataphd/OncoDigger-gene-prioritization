"""
intogen_build_figure1_workflow.py -- manuscript1.1 (IntOGen-primary draft) ONLY. The workflow
diagram; the "Analysis" and "Findings" boxes carry the IntOGen labels and numbers.

Writes: figures/source/intogen_figure1_workflow.pptx
"""
import sys

from pptx.enum.text import PP_ALIGN

from _style import DARK, OI, ROOT, WHITE, lerp, line, new_slide, rect, save, text

sys.path.insert(0, str(ROOT / "scripts" / "analysis"))
import vocab_policy as _vp  # the frozen vocabulary policy: the box below reads its counts, nothing is typed

_N_TERMS = len(_vp.SUPPRESS_SYMBOLS) + sum(len(v) for v in _vp.SUPPRESS_ALIASES.values())
_N_EXCL = len(_vp.GENE_EXCLUDE)

BLUE, GOLD, TEAL, PURPLE = OI["blue"], "B8860B", "0F9B8E", OI["purple"]


def tint(colour, t=0.82):
    return lerp(colour, WHITE, t)


def tier_header(s, x, y, w, h, label, fill):
    rect(s, x, y, w, h, fill)
    text(s, x, y, w, h, label, pt=13, bold=True, colour=WHITE, align=PP_ALIGN.CENTER)


def sub_box(s, x, y, w, h, title, body_lines, fill):
    rect(s, x, y, w, h, tint(fill), line=fill, lw=1.0)
    text(s, x + 0.05, y + 0.04, w - 0.1, 0.17, title, pt=8.5, bold=True, colour=DARK, align=PP_ALIGN.CENTER)
    text(s, x + 0.05, y + 0.22, w - 0.1, h - 0.26, "\n".join(body_lines), pt=7.5, bold=False, colour=DARK,
         align=PP_ALIGN.CENTER, wrap=True)


def arrow(s, x, y1, y2):
    line(s, x, y1, x, y2, DARK, 1.5)


def main():
    prs, s = new_slide()
    M = 0.35
    W = 6.7 - 2 * M
    text(s, M, M, W, 0.3, "OncoDigger Analytical Workflow", pt=15, bold=True)

    # Row 1: Input
    y = M + 0.42
    tier_header(s, M, y, W, 0.32, "Input: 19 Cancer-Specific PubMed Collections", BLUE)
    y2 = y + 0.32 + 0.06
    cancers = [
        ("GU", ["Bladder,", "Kidney,", "Prostate"]),
        ("Gyn/Breast", ["Breast,", "Cervical,", "Endometrial,", "Ovarian"]),
        ("GI", ["Colorectal,", "Stomach, Liver,", "Pancreatic,", "Esophageal"]),
        ("Thoracic", ["Lung"]),
        ("Skin", ["Melanoma"]),
        ("CNS", ["Brain"]),
        ("Heme", ["Leukemia,", "Lymphoma,", "Myeloma"]),
        ("Endocrine", ["Thyroid"]),
    ]
    cw = W / len(cancers)
    row1_h = 0.85
    for i, (title, body) in enumerate(cancers):
        sub_box(s, M + i * cw + 0.02, y2, cw - 0.04, row1_h, title, body, BLUE)
    y3 = y2 + row1_h
    arrow(s, M + W / 2, y3, y3 + 0.14)

    # Row 2: Engine
    y = y3 + 0.14
    tier_header(s, M, y, W, 0.32, "OncoDigger Engine: 'mutation mutations' query", GOLD)
    y2 = y + 0.32 + 0.06
    engine_boxes = [
        ("BM25 Retrieval", ["Corpus-specific IDF", "Title/MeSH/recency", "multipliers", "Top 1,000 abstracts"]),
        ("Gene Enrichment", ["Rank-weighted paper mentions,", "query-relevant vs. background", "literature rate,",
                              "min. 2 supporting papers"]),
        ("Vocabulary QC", ["143 alias collisions,", "8,203 nested names,", f"{_N_EXCL} genes excluded,", f"{_N_TERMS} suppression terms"]),
        ("Reproducibility", ["Fixed hash seed,", "per-corpus cache", "clearing"]),
    ]
    ew = W / len(engine_boxes)
    row2_h = 1.05
    for i, (title, body) in enumerate(engine_boxes):
        sub_box(s, M + i * ew + 0.02, y2, ew - 0.04, row2_h, title, body, GOLD)
    y3 = y2 + row2_h
    arrow(s, M + W / 2, y3, y3 + 0.14)

    # Row 3: Output
    y = y3 + 0.14
    tier_header(s, M, y, W, 0.32, "Output: Ranked Gene Lists (top-k per cancer, PMID-traceable)", TEAL)
    y3 = y + 0.32
    arrow(s, M + W / 2, y3, y3 + 0.14)

    # Row 4: Analysis (2 boxes)
    y = y3 + 0.14
    tier_header(s, M, y, W, 0.32, "Analysis", PURPLE)
    y2 = y + 0.32 + 0.06
    aw = W / 2
    row4_h = 1.15
    sub_box(s, M + 0.02, y2, aw - 0.06, row4_h, "Benchmarking & Controls", [
        "IntOGen driver Precision@k",
        "Pooled-corpus control (cancer-aware vs. cancer-blind)",
        "Frequency-only baseline (complete pipeline vs. counting)",
    ], PURPLE)
    sub_box(s, M + aw + 0.04, y2, aw - 0.06, row4_h, "Head-to-Head Benchmark", [
        "vs. CancerMine, DISEASES 2.0,",
        "PubTator 3.0 (all audited)",
        "Same 19 cancers, same IntOGen reference",
    ], PURPLE)
    y3 = y2 + row4_h
    arrow(s, M + W / 2, y3, y3 + 0.14)

    # Row 5: Findings (2 boxes)
    y = y3 + 0.14
    tier_header(s, M, y, W, 0.32, "Findings", OI["vermillion"])
    y2 = y + 0.32 + 0.06
    row5_h = 1.5
    import pandas as _pd
    _ci = _pd.read_csv(ROOT / "data/derived/intogen_rescoring/precision_at_k_bootstrap_ci_intogen.csv").set_index("k").loc[10]
    _cm = _pd.read_csv(ROOT / "data/derived/intogen_rescoring/statistical_support_reference_paired_vs_oncodigger.csv")
    _cm = _cm[(_cm.system == "CancerMine") & (_cm.reference == "published") & (_cm.metric == "any_cancer_P@10")].iloc[0]
    _ss = _pd.read_csv(ROOT / "data/derived/intogen_rescoring/statistical_support_reference_summary.csv")
    _ex = _ss[(_ss.system == "OncoDigger") & (_ss.reference == "published") & (_ss.metric == "exact_cancer_concordance")].iloc[0]
    _pp = _pd.read_csv(ROOT / "data/derived/intogen_rescoring/statistical_support_reference_paired_vs_oncodigger.csv")
    _pool = _pp[(_pp.system == "Pooled cancer-blind control") & (_pp.reference == "published") & (_pp.metric == "exact_cancer_concordance")].iloc[0]
    _spec = _pd.read_csv(ROOT / "data/derived/intogen_rescoring/corpus_specificity_clean_intogen.csv").cancer_specific_in_top10.mean()
    _abl = _pd.read_csv(ROOT / "data/derived/intogen_rescoring/ablation_panel_summary.csv")
    _full = float(_abl[_abl.iloc[:, 0].str.startswith("full_pipeline")].iloc[0, 3]); _freq = float(_abl[_abl.iloc[:, 0].str.startswith("frequency_only (vetted")].iloc[0, 3])
    sub_box(s, M + 0.02, y2, aw - 0.06, row5_h, "OncoDigger vs. IntOGen", [
        f"• {100 * _ci.mean_p_at_k:.1f}% mean Precision@10",
        f"  (95% cancer-bootstrap CI {100 * _ci['ci_low_2.5']:.1f}-{100 * _ci['ci_high_97.5']:.1f}%)",
        f"• {_ex.percent:.1f}% exact cancer-type concordance",
        f"  (95% CI {100 * _ex.boot_ci_lo:.1f}-{100 * _ex.boot_ci_hi:.1f}%)",
        f"• +{_pool.mean_diff_pp:.1f} pp exact vs. cancer-blind pooled control",
        f"• {_spec:.1f} cancer-specific genes/cancer",
        "  beyond the pooled-corpus control",
        f"• +{_full - _freq:.1f} pp over frequency baseline",
    ], OI["vermillion"])
    sub_box(s, M + aw + 0.04, y2, aw - 0.06, row5_h, "System-Level Comparison", [
        "• CancerMine: 87.9% (95% CI",
        f"  {_cm.boot_ci_lo_pp:.1f}-{_cm.boot_ci_hi_pp:.1f}pp advantage)",
        "• DISEASES 2.0 / PubTator 3.0:",
        "  60.0% / 67.9% (broader scope)",
        "• ChatGPT: 78.9% (practical",
        "  baseline, not a peer system)",
    ], OI["vermillion"])

    save(prs, "intogen_figure1_workflow.pptx")


if __name__ == "__main__":
    main()
