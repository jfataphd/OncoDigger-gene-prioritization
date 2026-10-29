"""
intogen_build_figure5_literature_weight.py -- manuscript1.1 (IntOGen-primary draft) ONLY.
This is manuscript1.1's Figure 5, "literature weight of unrecovered driver annotations", scored
against IntOGen's 1,740 gene-cancer driver pairs.

Inputs: data/derived/intogen_rescoring/intogen_literature_weight_audit.csv
Writes: figures/source/intogen_figure5_literature_weight.pptx
"""
import math

import pandas as pd
from pptx.enum.text import PP_ALIGN

from _style import DARK, GREY, LIGHT, M, OI, ROOT, axes, line, marker, new_slide, panel, rect, save, text

BAND_COLOURS = {"top10": OI["blue"], "ranks11_50": OI["orange"], "beyond50": OI["vermillion"]}
BAND_LABELS = {"top10": "Top 10", "ranks11_50": "Ranks 11–50", "beyond50": "Ranked beyond 50"}


def band_of(r):
    if pd.isna(r):
        return "not_surfaced"
    if r <= 10:
        return "top10"
    if r <= 50:
        return "ranks11_50"
    return "beyond50"


def main():
    d = pd.read_csv(ROOT / "data" / "derived" / "intogen_rescoring" / "intogen_literature_weight_audit.csv")
    d["band"] = d["engine_rank"].apply(band_of)
    total = len(d)
    from scipy import stats as _st
    _s = d[d["engine_rank"].notna()]
    _ns = len(_s)
    _rho = _st.spearmanr(_s["engine_rank"], _s["supporting_papers"])[0]
    _ps = pd.read_csv(ROOT / "data/derived/intogen_rescoring/intogen_full_corpus_presence_summary.csv").set_index("category")["n"]
    _np, _na = int(_ps["present_but_outside_window"]), int(_ps["truly_absent_from_corpus"])

    prs, s = new_slide()
    text(s, M, M, 6.0, 0.25, "Literature support and recovery of IntOGen driver annotations", pt=11, bold=True)

    # (A) box summary by band, log10 y
    panel(s, M, M + 0.35, "A", f"Supporting papers by recovery band ({int(d['engine_rank'].notna().sum()):,} of {total:,} IntOGen driver pairs surfaced)")
    x0, y0, w, h = M + 1.0, M + 0.85, 4.0, 2.6
    yticks = [0, 1, 2, 3]
    P = axes(s, x0, y0, w, h, (0, 4), (0, 3), [0.5, 1.5, 2.5, 3.5], yticks,
             "", "Supporting papers (log scale)", fmt_x=lambda v: "", fmt_y=lambda v: f"{10**v:,.0f}", grid=True)

    order = ["top10", "ranks11_50", "beyond50"]
    bw = 0.7
    for i, band in enumerate(order):
        sub = d[d["band"] == band]["supporting_papers"]
        q25, med, q75 = sub.quantile(.25), sub.median(), sub.max()
        lo, hi = sub.min(), sub.max()
        cx = 0.5 + i + 1.0
        col = BAND_COLOURS[band]
        l25, l75 = math.log10(sub.quantile(.25)), math.log10(sub.quantile(.75))
        lmed = math.log10(sub.median())
        llo, lhi = math.log10(sub.min()), math.log10(sub.max())
        xpix = P(cx, 0)[0]
        # whiskers
        line(s, xpix, P(cx, llo)[1], xpix, P(cx, l25)[1], DARK, 1.0)
        line(s, xpix, P(cx, l75)[1], xpix, P(cx, lhi)[1], DARK, 1.0)
        # box
        bx0 = P(cx - bw / 2, 0)[0]
        by0 = P(0, l75)[1]
        bh = P(0, l25)[1] - P(0, l75)[1]
        rect(s, bx0, by0, bw, bh, col, line=DARK, lw=0.75)
        # median line
        mx0, my = P(cx - bw / 2, lmed)
        line(s, mx0, my, mx0 + bw, my, "FFFFFF", 1.75)
        text(s, xpix - 0.4, P(cx, 0)[1] + 0.06, 0.8, 0.16, BAND_LABELS[band], pt=7.5, align=PP_ALIGN.CENTER)
        n = len(sub)
        text(s, xpix - 0.4, P(cx, llo)[1] + 0.06, 0.8, 0.14, f"n={n}", pt=7, colour=GREY, align=PP_ALIGN.CENTER)

    # never-surfaced callout
    nsurf = (d["band"] == "not_surfaced").sum()
    cbw, cbh = 1.7, 0.92
    cbx = x0 + w - cbw - 0.05
    cby0 = y0 - 0.12  # nudged up from y0 + 0.05 so the box clears the top10 whisker below it
    rect(s, cbx, cby0, cbw, cbh, "F3F3F3", line=GREY, lw=0.5)
    text(s, cbx + 0.05, cby0 + 0.05, cbw - 0.1, 0.18, f"{nsurf}/{total} ({100*nsurf/total:.1f}%)", pt=9, bold=True,
         align=PP_ALIGN.CENTER)
    text(s, cbx + 0.05, cby0 + 0.27, cbw - 0.1, 0.18, "never surfaced in window\n(no value to plot)", pt=7,
         colour=GREY, align=PP_ALIGN.CENTER)
    text(s, cbx + 0.05, cby0 + 0.57, cbw - 0.1, 0.30, f"of which {_np} present\nelsewhere, {_na} no lexicon-\ndetected mention (Suppl. Fig. S3)",
         pt=6.5, colour=GREY, align=PP_ALIGN.CENTER)

    # (B) rank vs support scatter, log-log
    top = M + 3.75
    n_surfaced = d["engine_rank"].notna().sum()
    panel(s, M, top, "B", f"OncoDigger rank vs. supporting papers ({n_surfaced} surfaced IntOGen pairs)")
    x0, y0, w, h = M + 1.0, top + 0.45, 3.85, 2.9
    surfaced = d.dropna(subset=["engine_rank"]).copy()
    surfaced["log_rank"] = surfaced["engine_rank"].apply(math.log10)
    surfaced["log_sp"] = surfaced["supporting_papers"].apply(lambda v: math.log10(max(v, 1)))
    P2 = axes(s, x0, y0, w, h, (0, 2.6), (0, 3), [0, 1, 2], [0, 1, 2, 3],
              "OncoDigger rank (log scale)", "Supporting papers (log scale)",
              fmt_x=lambda v: f"{10**v:,.0f}", fmt_y=lambda v: f"{10**v:,.0f}")
    for _, r in surfaced.iterrows():
        col = BAND_COLOURS.get(r["band"], GREY)
        px, py = P2(r["log_rank"], r["log_sp"])
        marker(s, px, py, 0.05, col)

    # correlation callout
    cbw, cbh = 2.1, 0.62
    cbx, cby = x0 + w - cbw - 0.05, y0 + 0.05
    rect(s, cbx, cby, cbw, cbh, "F3F8FC", line=OI["sky"], lw=0.75)
    text(s, cbx + 0.06, cby + 0.05, cbw - 0.1, 0.15, f"Spearman ρ = {_rho:.3f} (n={_ns})".replace("-", "−"), pt=7.5, bold=True)
    text(s, cbx + 0.06, cby + 0.24, cbw - 0.1, 0.15, "Same ordering (top10 > 11–50 > beyond 50)", pt=7)
    text(s, cbx + 0.06, cby + 0.41, cbw - 0.1, 0.15, "held in all 19/19 study cancers", pt=7)

    # legend
    legend_items = order
    ly0 = y0 + h + 0.48
    for i, key in enumerate(legend_items):
        xx = x0 + i * 1.3
        yy = ly0
        marker(s, xx + 0.07, yy + 0.08, 0.08, BAND_COLOURS[key])
        text(s, xx + 0.2, yy, 1.1, 0.16, BAND_LABELS[key], pt=7)

    save(prs, "intogen_figure5_literature_weight.pptx")


if __name__ == "__main__":
    main()
