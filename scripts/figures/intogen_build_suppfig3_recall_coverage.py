"""
intogen_build_suppfig3_recall_coverage.py -- manuscript1.1 (IntOGen-primary draft) ONLY. Supplementary Figure S3: why most IntOGen driver pairs
do not reach the ranking, in three panels. It replaces and combines the former Supplementary Figures S6 and S7 (retired 2026-10-18; their
scripts and images are not part of the final supplement).

(A) All 1,740 IntOGen driver pairs split into surfaced in the standard analysis window, present in the full corpus but outside the window, and
    no lexicon-detected mention anywhere in the full corpus (Supplementary Results S10).
(B) Full-corpus mention count of the present-but-outside-window pairs (excludes the confirmed false-alias collision WAS/Brain Cancer), log-binned.
(C) Recall at ranks 10, 20, 50 and any rank under three progressively stricter filters (Supplementary Results S11).

Inputs: data/derived/intogen_rescoring/intogen_full_corpus_presence_check.csv, intogen_full_corpus_presence_summary.csv,
        intogen_cleaned_recall_summary.csv
Writes: figures/source/intogen_suppfig3_recall_coverage.pptx
"""
import pandas as pd
from pptx.enum.text import PP_ALIGN

from _style import DARK, GREY, LIGHT, M, OI, ROOT, line, panel, rect, save, new_slide, text

SURF, PRESENT, ABSENT = OI["green"], OI["sky"], OI["purple"]
FCOL = [OI["grey"], OI["orange"], OI["blue"]]
FLAB = ["No filter", "Remove zero-mention pairs", "Also require 2 or more IntOGen methods"]
DEPTHS = [("rank10_pct", "Rank 10"), ("rank20_pct", "Rank 20"), ("rank50_pct", "Rank 50"), ("any_rank_pct", "Any rank")]


def main():
    der = ROOT / "data" / "derived" / "intogen_rescoring"
    det = pd.read_csv(der / "intogen_full_corpus_presence_check.csv")
    summ = pd.read_csv(der / "intogen_full_corpus_presence_summary.csv").set_index("category")
    rec = pd.read_csv(der / "intogen_cleaned_recall_summary.csv")
    n_surf = int(summ.loc["surfaced_in_window", "n"])
    n_present = int(summ.loc["present_but_outside_window", "n"])
    n_absent = int(summ.loc["truly_absent_from_corpus", "n"])
    total = n_surf + n_present + n_absent

    prs, s = new_slide()
    text(s, M, M, 6.0, 0.25, "Why most IntOGen driver pairs do not reach the ranking", pt=11, bold=True)

    # ---- (A) stacked bar ----
    panel(s, M, M + 0.35, "A", f"All {total:,} IntOGen driver pairs")
    x0, y0, w, bh = M + 0.5, M + 0.95, 5.6, 0.45
    segs = [("Surfaced", n_surf, SURF), ("Present, outside window", n_present, PRESENT),
            ("No detected mention", n_absent, ABSENT)]
    x = x0
    for lab, n, col in segs:
        sw = w * n / total
        rect(s, x, y0, sw, bh, col)
        if sw > 0.7:
            text(s, x, y0 + bh / 2 - 0.1, sw, 0.2, f"{n} ({100 * n / total:.1f}%)", pt=8.5, bold=True, colour="FFFFFF",
                 align=PP_ALIGN.CENTER)
        x += sw
    ly = y0 + bh + 0.12
    lx = x0
    for lab, n, col in segs:
        rect(s, lx, ly, 0.14, 0.14, col)
        text(s, lx + 0.2, ly - 0.03, 2.2, 0.2, f"{lab}: {n} ({100 * n / total:.1f}%)", pt=7.5)
        lx += 1.95
    text(s, M, ly + 0.25, 6.0, 0.4,
         f"\"Never surfaced\" ({n_present + n_absent}) is not \"absent from the literature\": {n_absent} of them "
         f"({100 * n_absent / (n_present + n_absent):.1f}%) have no lexicon-detected mention, which reflects vocabulary detection, not confirmed absence.",
         pt=7.5, italic=True, colour=GREY, wrap=True)

    # ---- (B) histogram ----
    top = M + 2.55
    present_df = det[(det["full_corpus_postings"] > 0) & (~det["false_alias_excluded"])]
    panel(s, M, top, "B", f"Full-corpus mentions, present-but-outside-window pairs (n={len(present_df)})")
    bins = [(1, 2, "1"), (2, 5, "2-4"), (5, 10, "5-9"), (10, 20, "10-19"), (20, 50, "20-49"), (50, 100, "50-99"), (100, 10**9, "100+")]
    counts = [(lab, int(((present_df["full_corpus_postings"] >= lo) & (present_df["full_corpus_postings"] < hi)).sum())) for lo, hi, lab in bins]
    bx0, by0, bw, bh2 = M + 0.6, top + 0.55, 5.4, 1.45
    maxn = max(n for _, n in counts)
    cw = bw / len(counts)
    for t in (0, 50, 100, 150, 200):
        if t > maxn * 1.15:
            continue
        yy = by0 + bh2 - t / (maxn * 1.15) * bh2
        line(s, bx0, yy, bx0 + bw, yy, LIGHT, 0.5)
        text(s, bx0 - 0.4, yy - 0.08, 0.35, 0.16, str(t), pt=7.5, align=PP_ALIGN.RIGHT)
    line(s, bx0, by0, bx0, by0 + bh2, DARK, 0.75)
    line(s, bx0, by0 + bh2, bx0 + bw, by0 + bh2, DARK, 0.75)
    for i, (lab, n) in enumerate(counts):
        h_i = n / (maxn * 1.15) * bh2
        xx = bx0 + i * cw
        rect(s, xx + cw * 0.12, by0 + bh2 - h_i, cw * 0.76, h_i, PRESENT)
        text(s, xx, by0 + bh2 - h_i - 0.17, cw, 0.16, str(n), pt=7.5, bold=True, align=PP_ALIGN.CENTER)
        text(s, xx, by0 + bh2 + 0.03, cw, 0.16, lab, pt=7.5, align=PP_ALIGN.CENTER)
    med = int(present_df["full_corpus_postings"].median())
    text(s, bx0, by0 + bh2 + 0.22, bw, 0.16, f"Full-corpus mentions (median = {med})", pt=8.5, align=PP_ALIGN.CENTER)
    text(s, bx0 - 0.4 - 0.1 - bh2 / 2, by0 + bh2 / 2 - 0.09, bh2, 0.18, "Number of pairs", pt=8.5, align=PP_ALIGN.CENTER, rotate=270)

    # ---- (C) recall under filters ----
    topc = by0 + bh2 + 0.55
    panel(s, M, topc, "C", "Recall under progressively stricter, independently justified filters")
    cx0, cy0, cw_, ch = M + 0.85, topc + 0.5, 5.2, 2.15
    for t in (0, 20, 40, 60, 80):
        yy = cy0 + ch - t / 80 * ch
        line(s, cx0, yy, cx0 + cw_, yy, LIGHT, 0.5)
        text(s, cx0 - 0.55, yy - 0.08, 0.5, 0.16, f"{t}%", pt=7.5, align=PP_ALIGN.RIGHT)
    line(s, cx0, cy0, cx0, cy0 + ch, DARK, 0.75)
    line(s, cx0, cy0 + ch, cx0 + cw_, cy0 + ch, DARK, 0.75)
    gw = cw_ / len(DEPTHS)
    bwd = gw * 0.26
    for gi, (col, lab) in enumerate(DEPTHS):
        gx = cx0 + gi * gw
        for fi, frow in rec.iterrows():
            val = frow[col]
            bhh = val / 80 * ch
            xx = gx + gw / 2 - 1.5 * bwd + fi * bwd
            rect(s, xx, cy0 + ch - bhh, bwd * 0.88, bhh, FCOL[fi])
            text(s, xx - 0.05, cy0 + ch - bhh - 0.15, bwd + 0.1, 0.14, f"{val:.0f}%", pt=6, bold=True, colour=FCOL[fi], align=PP_ALIGN.CENTER)
        text(s, gx, cy0 + ch + 0.04, gw, 0.18, lab, pt=8.5, align=PP_ALIGN.CENTER)
    text(s, cx0 - 0.55 - 0.1 - ch / 2, cy0 + ch / 2 - 0.09, ch, 0.18, "Recall", pt=9, align=PP_ALIGN.CENTER, rotate=270)
    ly = cy0 + ch + 0.32
    lx = cx0 - 0.5
    for col, lab, wl in zip(FCOL, FLAB, (0.95, 1.9, 2.4)):
        rect(s, lx, ly, 0.14, 0.14, col)
        text(s, lx + 0.2, ly - 0.03, wl, 0.2, lab, pt=7.5)
        lx += wl + 0.25
    save(prs, "intogen_suppfig3_recall_coverage.pptx")


if __name__ == "__main__":
    main()
