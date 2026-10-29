"""
intogen_build_figure4_head_to_head.py -- manuscript1.1 (IntOGen-primary draft) ONLY.
Scores OncoDigger, CancerMine, DISEASES and PubTator against IntOGen's any-cancer driver
set (this is manuscript1.1's Figure 4, "comparison with
published literature-mining resources").

Inputs: data/derived/intogen_rescoring/oncodigger_intogen_summary.csv,
        cancermine_intogen_summary.csv, diseases_intogen_summary.csv,
        pubtator_v2_intogen_summary.csv, figure4_headtohead_pvalues.csv,
        pubtator_original_intogen_detail.csv
Writes: figures/source/intogen_figure4_head_to_head.pptx
"""
import pandas as pd
from pptx.enum.text import PP_ALIGN

from _style import DARK, GREY, LIGHT, M, OI, ROOT, W, line, marker, new_slide, panel, rect, save, text

OD = OI["blue"]
CM = OI["vermillion"]
DS = OI["green"]
PT = OI["purple"]
UNAUDITED = OI["grey"]


def load_summary(name):
    der = ROOT / "data" / "derived" / "intogen_rescoring"
    s = pd.read_csv(der / f"{name}_intogen_summary.csv").set_index("cancer")
    col = "p10_intogen_any" if "p10_intogen_any" in s.columns else "intogen_any"
    return s[col]


def main():
    der = ROOT / "data" / "derived" / "intogen_rescoring"
    od = load_summary("oncodigger")
    cm = load_summary("cancermine")
    ds = load_summary("diseases")
    pt = load_summary("pubtator_v2")
    # PubTator's pre-correction percentage is not shown: the source file it was computed
    # from has no traceable generator script, and a live reconstruction of the original
    # single-parent-MeSH methodology (pubtator_v1_single_parent_reconstruction.py) returned
    # 67.9% -- indistinguishable from the corrected value, not the ~65.8% once reported --
    # confirming PubTator's live relations data has changed since the original benchmark
    # window and the original figure cannot be verified (Supplementary Table S3).

    prs, s = new_slide()

    # (A) Mean Precision@10 grouped bars
    panel(s, M, M, "A", "Mean Precision@10 vs IntOGen (any cancer type)")
    x0, y0, w, h = M + 0.85, M + 0.68, 5.0, 2.1
    for t in (0, 0.25, 0.5, 0.75, 1.0):
        yy = y0 + h - t * h
        line(s, x0, yy, x0 + w, yy, LIGHT, 0.5)
        text(s, x0 - 0.55, yy - 0.08, 0.5, 0.16, f"{t:.0%}", pt=8, align=PP_ALIGN.RIGHT)
    line(s, x0, y0, x0, y0 + h, DARK, 0.75)
    line(s, x0, y0 + h, x0 + w, y0 + h, DARK, 0.75)

    # CI label is the 95% paired cancer-bootstrap CI on the OncoDigger-minus-tool
    # difference (Supplementary Results S4) -- treated as the primary inferential evidence,
    # with the paired t-test p-value shown smaller, beneath it, as secondary corroboration
    # (matching the paired-bootstrap-primary framing already used in Methods/Results text).
    # intervals are read from the canonical paired-bootstrap output (never typed in)
    _pp = pd.read_csv(ROOT / "data" / "derived" / "intogen_rescoring" / "statistical_support_reference_paired_vs_oncodigger.csv")
    _pp = _pp[(_pp.reference == "published") & (_pp.metric == "any_cancer_P@10")].set_index("system")

    def _ci(sysname):
        r = _pp.loc[sysname]
        return f"Δ {r.boot_ci_lo_pp:.1f} to {r.boot_ci_hi_pp:.1f} pp".replace("Δ -0.0", "Δ −0.0")

    # paired t-test p-values come from figure4_headtohead_pvalues.csv (intogen_figure4_headtohead_data.py), never typed in
    _pv = pd.read_csv(ROOT / "data" / "derived" / "intogen_rescoring" / "figure4_headtohead_pvalues.csv").set_index("tool").p_ttest

    def _p(tool):
        v = float(_pv.loc[tool])
        return "p < 0.0001" if v < 1e-4 else f"p = {v:.3f}"

    bars = [
        ("OncoDigger", od.mean(), OD, "(19 cancers)", None, None),
        ("CancerMine", cm.mean(), CM, "(19; audited)", _ci("CancerMine"), _p("cancermine")),
        ("DISEASES\n2.0", ds.mean(), DS, "(19)", _ci("DISEASES 2.0"), _p("diseases")),
        ("PubTator\n3.0", pt.mean(), PT, "(19; subtype-aggregated)", _ci("PubTator 3.0"), _p("pubtator_v2")),
    ]
    bw = 0.62
    gap = (w - bw * len(bars)) / (len(bars) + 1)
    bx = x0 + gap
    for label, val, col, n_lab, ci_lab, p_lab in bars:
        bh = val * h
        rect(s, bx, y0 + h - bh, bw, bh, col)
        if val >= 0.85:
            text(s, bx, y0 + h - bh + 0.06, bw, 0.16, f"{val:.1%}", pt=8, bold=True, colour="FFFFFF",
                 align=PP_ALIGN.CENTER)
            top = y0 + h - bh
        else:
            text(s, bx - 0.1, y0 + h - bh - 0.19, bw + 0.2, 0.16, f"{val:.1%}", pt=8, bold=True, colour=col,
                 align=PP_ALIGN.CENTER)
            top = y0 + h - bh - 0.19
        if ci_lab:
            text(s, bx - 0.25, top - 0.18, bw + 0.5, 0.16, ci_lab, pt=7.5, bold=True, colour=DARK,
                 align=PP_ALIGN.CENTER)
            text(s, bx - 0.25, top - 0.33, bw + 0.5, 0.15, p_lab, pt=8, italic=True, colour="444444",
                 align=PP_ALIGN.CENTER)
        text(s, bx - 0.15, y0 + h + 0.04, bw + 0.3, 0.32, label, pt=8, bold=True, colour=col, align=PP_ALIGN.CENTER)
        text(s, bx - 0.15, y0 + h + 0.34, bw + 0.3, 0.16, n_lab, pt=7.0, align=PP_ALIGN.CENTER)
        bx += bw + gap
    text(s, x0 - 0.55 - 0.1 - h / 2, y0 + h / 2 - 0.09, h, 0.18, "Mean P@10", pt=10, align=PP_ALIGN.CENTER,
         rotate=270)
    text(s, x0, y0 + h + 0.58, w, 0.16,
         "Δ: 95% paired cancer-bootstrap CI vs OncoDigger (primary); p-values: paired t-test (secondary)",
         pt=7, italic=True, colour=GREY, align=PP_ALIGN.CENTER)

    # (B) Per-cancer dumbbell, OncoDigger vs CancerMine
    top = M + 3.55
    panel(s, M, top, "B", "Per-cancer Precision@10: OncoDigger vs CancerMine (audited)")
    df = pd.DataFrame({"OD": od, "CM": cm}).dropna()
    df["d"] = df["OD"] - df["CM"]
    df = df.sort_values(["d"], ascending=False)
    bx, by, bwid = M + 1.45, top + 0.55, 3.9
    rh = 0.13
    bh_tot = rh * len(df)
    for t in (0, 0.25, 0.5, 0.75, 1.0):
        xx = bx + t * bwid
        line(s, xx, by, xx, by + bh_tot, LIGHT, 0.5)
        text(s, xx - 0.25, by + bh_tot + 0.04, 0.5, 0.16, f"{t:.0%}", pt=8, align=PP_ALIGN.CENTER)
    line(s, bx, by + bh_tot, bx + bwid, by + bh_tot, DARK, 0.75)
    line(s, bx, by, bx, by + bh_tot, DARK, 0.75)
    text(s, bx, by + bh_tot + 0.22, bwid, 0.18, "Precision@10", pt=10, align=PP_ALIGN.CENTER)
    for i, (cancer, r) in enumerate(df.iterrows()):
        yy = by + (i + 0.5) * rh
        text(s, M, yy - 0.08, 1.4, 0.16, cancer, pt=7, align=PP_ALIGN.RIGHT)
        xo, xc = bx + r["OD"] * bwid, bx + r["CM"] * bwid
        line(s, xc, yy, xo, yy, GREY, 1.0)
        marker(s, xc, yy, 0.075, CM)
        marker(s, xo, yy, 0.075, OD)
        d_pp = r["d"] * 100
        sign = "+" if d_pp >= 0 else "−"
        text(s, bx + bwid + 0.08, yy - 0.08, 0.5, 0.16, f"{sign}{abs(d_pp):.0f}", pt=7, bold=True,
             colour=(OD if d_pp >= 0 else CM))
    text(s, bx + bwid + 0.02, by - 0.2, 0.7, 0.16, "Δ (pp)", pt=8, bold=True)
    for i, (lab, col) in enumerate((("OncoDigger", OD), ("CancerMine", CM))):
        marker(s, bx + 0.15 + i * 1.3, top + 0.36, 0.075, col)
        text(s, bx + 0.25 + i * 1.3, top + 0.28, 1.1, 0.16, lab, pt=8, bold=True, colour=col)

    # (C) Auditing changes the answer
    top = M + 7.0
    panel(s, M, top, "C", "Effect of auditing: CancerMine before and after; PubTator final")
    cx, cy, cw, ch = M + 1.3, top + 0.35, 3.5, 0.55
    for t in (0, 0.25, 0.5, 0.75, 1.0):
        xx = cx + t * cw
        line(s, xx, cy, xx, cy + ch, LIGHT, 0.5)
        text(s, xx - 0.25, cy + ch + 0.05, 0.5, 0.16, f"{t:.0%}", pt=7, align=PP_ALIGN.CENTER)
    line(s, cx, cy + ch, cx + cw, cy + ch, DARK, 0.75)
    line(s, cx, cy, cx, cy + ch, DARK, 0.75)

    rows = [
        ("CancerMine", 0.842, cm.mean(), CM,
         "PTGDR artifact (1,329 mis-tagged “-AS1” sentences) + DOID mapping errors incl. Brain Cancer subtypes"),
        ("PubTator 3.0", None, pt.mean(), PT,
         "Single-parent-MeSH mapping missed subtype-tagged relations; fixing it also resolved Liver/ALS contamination"),
    ]
    rh2 = ch / 2
    for i, (label, before, after, col, note) in enumerate(rows):
        yy = cy + i * rh2
        bar_h = rh2 * 0.34
        if before is not None:
            rect(s, cx, yy + bar_h * 0.15, before * cw, bar_h, UNAUDITED)
            text(s, cx + before * cw + 0.05, yy + bar_h * 0.15 - 0.02, 0.5, bar_h + 0.04, f"{before:.1%}", pt=7,
                 colour=GREY)
        rect(s, cx, yy + bar_h * 1.55, after * cw, bar_h, col)
        text(s, cx + after * cw + 0.05, yy + bar_h * 1.55 - 0.02, 0.5, bar_h + 0.04, f"{after:.1%}", pt=7, bold=True,
             colour=col)
        text(s, M, yy - 0.02, 1.25, rh2, label, pt=7.5, bold=True, align=PP_ALIGN.RIGHT, colour=col)
    text(s, cx, cy + ch + 0.22, cw, 0.16, "Mean Precision@10", pt=9, align=PP_ALIGN.CENTER)
    ly = cy + ch + 0.42
    rect(s, cx + 0.3, ly + 0.02, 0.11, 0.11, UNAUDITED)
    text(s, cx + 0.46, ly - 0.01, 1.1, 0.16, "Unaudited pass", pt=7)
    rect(s, cx + 1.75, ly + 0.02, 0.09, 0.11, CM)
    rect(s, cx + 1.84, ly + 0.02, 0.09, 0.11, PT)
    text(s, cx + 1.98, ly - 0.01, 1.5, 0.16, "Audited (this study)", pt=7)
    note_y = cy + ch + 0.60
    note_w = W - 2 * M
    for label, before, after, col, note in rows:
        text(s, M, note_y, note_w, 0.3, f"• {note}", pt=8, italic=True, colour=DARK, wrap=True)
        lines = -(-len(note) // 80)  # ceil div: ~80 chars/line at this width/font
        note_y += 0.15 * lines

    save(prs, "intogen_figure4_head_to_head.pptx")


if __name__ == "__main__":
    main()
