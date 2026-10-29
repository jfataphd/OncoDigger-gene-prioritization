"""
intogen_build_figure3_validation.py -- manuscript1.1 (IntOGen-primary draft) ONLY.
All three panels are scored against IntOGen's any-cancer driver set.

Inputs: data/derived/intogen_rescoring/freq_baseline_per_cancer_detail_intogen.csv,
        freq_baseline_precision_at_k_aliases_intogen.csv, corpus_specificity_clean_intogen.csv
Writes: figures/source/intogen_figure3_validation.pptx
"""
import pandas as pd
from pptx.enum.dml import MSO_LINE_DASH_STYLE
from pptx.enum.text import PP_ALIGN

from _style import (DARK, GREY, LIGHT, M, OI, ROOT, axes, line, marker, new_slide, panel, polygon, polyline, rect,

                    save, text)

import pandas as _pd
# distinct IntOGen driver genes (633) over the 19,250 protein-coding genes (Seal et al., 2026): the random baseline
_N_DRIVERS = _pd.read_csv(ROOT / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv",
                          sep=chr(9))["SYMBOL"].nunique()

KS = [5, 10, 15, 20, 30, 40, 50]
OD, FQ = OI["blue"], OI["vermillion"]


def main():
    der = ROOT / "data" / "derived" / "intogen_rescoring"
    det = pd.read_csv(der / "freq_baseline_per_cancer_detail_intogen.csv")
    summ = pd.read_csv(der / "freq_baseline_precision_at_k_aliases_intogen.csv").set_index("k")
    spec = pd.read_csv(der / "corpus_specificity_clean_intogen.csv")
    prs, s = new_slide()

    # (A)
    panel(s, M, M, "A", "BM25 + enrichment vs frequency-only baseline: Precision@k")
    x0, y0, w, h = M + 0.75, M + 0.42, 5.1, 2.05
    P = axes(s, x0, y0, w, h, (5, 50), (0, 1), KS, [0, 0.25, 0.5, 0.75, 1.0], "Gene-rank cutoff (k)",
             "Mean P@k (± 1 SD)", fmt_y=lambda v: f"{v:.0%}")
    for pre, col in (("OD", OD), ("Freq", FQ)):
        mu = [det[f"{pre}@{k}"].mean() for k in KS]
        sd = [det[f"{pre}@{k}"].std() for k in KS]
        up = [P(k, min(1, m + e)) for k, m, e in zip(KS, mu, sd)]
        lo = [P(k, max(0, m - e)) for k, m, e in zip(KS, mu, sd)]
        polygon(s, up + lo[::-1], col, alpha=0.18)
        pts = [P(k, m) for k, m in zip(KS, mu)]
        polyline(s, pts, col, 1.5)
        for px, py in pts:
            marker(s, px, py, 0.08, col)
    for k in (10, 30):
        a, b = P(k, summ.loc[k, "OD_mean"]), P(k, summ.loc[k, "Freq_mean"])
        line(s, a[0], a[1], b[0], b[1], DARK, 0.75)
        text(s, a[0] + 0.05, (a[1] + b[1]) / 2 - 0.09, 1.1, 0.18, f"Δ = +{summ.loc[k, 'delta_pp']:.1f} pp", pt=8,
             bold=True)
    ry = P(5, _N_DRIVERS / 19250)[1]
    line(s, x0, ry, x0 + w, ry, DARK, 0.75, MSO_LINE_DASH_STYLE.DASH)
    text(s, x0 + 0.1, ry - 0.2, 1.5, 0.16, f"random = {_N_DRIVERS / 19250:.1%}", pt=8)
    lx, ly = x0 + w - 2.3, y0 + h - 0.72
    rect(s, lx, ly, 2.25, 0.5, "FFFFFF", line=GREY)
    for i, (lab, col) in enumerate((("OncoDigger (BM25 + enrichment)", OD), ("Frequency-only baseline", FQ))):
        yy = ly + 0.13 + i * 0.22
        line(s, lx + 0.08, yy, lx + 0.32, yy, col, 1.5)
        marker(s, lx + 0.2, yy, 0.07, col)
        text(s, lx + 0.38, yy - 0.08, 1.7, 0.16, lab, pt=8, bold=True, colour=col)

    # (B)
    top = M + 2.95
    panel(s, M, top, "B", "Per-cancer Precision@10: OncoDigger vs frequency-only")
    det["d10"] = det["Δ@10"]
    b = det.sort_values(["d10", "cancer"], ascending=[False, True]).reset_index(drop=True)
    bx, by, bw = M + 1.45, top + 0.55, 3.9
    rh = 0.115
    bh = rh * len(b)
    for t in (0, 0.25, 0.5, 0.75, 1.0):
        xx = bx + t * bw
        line(s, xx, by, xx, by + bh, LIGHT, 0.5)
        text(s, xx - 0.25, by + bh + 0.04, 0.5, 0.16, f"{t:.0%}", pt=8, align=PP_ALIGN.CENTER)
    line(s, bx, by + bh, bx + bw, by + bh, DARK, 0.75)
    line(s, bx, by, bx, by + bh, DARK, 0.75)
    text(s, bx, by + bh + 0.22, bw, 0.18, "Precision@10", pt=10, align=PP_ALIGN.CENTER)
    for i, r in b.iterrows():
        yy = by + (i + 0.5) * rh
        text(s, M, yy - 0.08, 1.4, 0.16, r["cancer"], pt=7, align=PP_ALIGN.RIGHT)
        xo, xf = bx + r["OD@10"] * bw, bx + r["Freq@10"] * bw
        line(s, xf, yy, xo, yy, GREY, 1.0)
        marker(s, xf, yy, 0.075, FQ)
        marker(s, xo, yy, 0.075, OD)
        text(s, bx + bw + 0.08, yy - 0.08, 0.45, 0.16, f"+{r['d10']:.0f}", pt=7, bold=True, colour=OD)
    text(s, bx + bw + 0.02, by - 0.2, 0.6, 0.16, "Δ (pp)", pt=8, bold=True)
    for i, (lab, col) in enumerate((("OncoDigger", OD), ("Frequency-only", FQ))):
        marker(s, bx + 0.15 + i * 1.1, top + 0.36, 0.075, col)
        text(s, bx + 0.25 + i * 1.1, top + 0.28, 1.0, 0.16, lab, pt=8, bold=True, colour=col)

    # (C)
    top = M + 6.05
    panel(s, M, top, "C", "Corpus separation recovers genes absent from the pooled ranking")
    c = spec.sort_values(["cancer_specific_in_top10", "cancer_specific_intogen", "cancer"],
                         ascending=[False, False, True]).reset_index(drop=True)
    cx, cy, cw = M + 1.45, top + 0.4, 3.9
    rh = 0.1
    ch = rh * len(c)
    for t in range(0, 11, 2):
        xx = cx + t / 10 * cw
        line(s, xx, cy, xx, cy + ch, LIGHT, 0.5)
        text(s, xx - 0.2, cy + ch + 0.04, 0.4, 0.16, str(t), pt=8, align=PP_ALIGN.CENTER)
    line(s, cx, cy + ch, cx + cw, cy + ch, DARK, 0.75)
    line(s, cx, cy, cx, cy + ch, DARK, 0.75)
    text(s, cx, cy + ch + 0.22, cw, 0.18, "OncoDigger top-10 genes absent from the pooled-corpus top-10", pt=10,
         align=PP_ALIGN.CENTER)
    for i, r in c.iterrows():
        yy = cy + i * rh
        t1, n = r["cancer_specific_intogen"], r["cancer_specific_in_top10"]
        text(s, M, yy - 0.01, 1.4, rh, r["cancer"], pt=7, align=PP_ALIGN.RIGHT)
        rect(s, cx, yy + 0.015, t1 / 10 * cw, rh - 0.03, OD)
        if n > t1:
            rect(s, cx + t1 / 10 * cw, yy + 0.015, (n - t1) / 10 * cw, rh - 0.03, OI["orange"])
        text(s, cx + n / 10 * cw + 0.05, yy - 0.01, 0.3, rh, str(n), pt=7, bold=True, colour=OD)
    lx = cx + cw - 1.0   # inside the plot, bottom right (bars there are short)
    for i, (lab, col) in enumerate((("IntOGen driver", OD), ("Not IntOGen driver", OI["orange"]))):
        rect(s, lx, cy + ch - 0.42 + i * 0.18, 0.11, 0.11, col)
        text(s, lx + 0.15, cy + ch - 0.45 + i * 0.18, 0.85, 0.16, lab, pt=7)
    save(prs, "intogen_figure3_validation.pptx")


if __name__ == "__main__":
    main()
