"""
intogen_build_suppfig2_confidence_breakdown.py -- manuscript1.1 (IntOGen-primary draft) ONLY.
Supplementary figure comparing IntOGen's own
internal confidence signals between the 754 driver pairs OncoDigger's ranking surfaced and the
986 it never surfaced (Results 5.5).

Single panel: distribution of IntOGen method-count agreement (1-7), surfaced vs never-surfaced
pairs. The median-q-value / median-%-cohort-mutated panel (formerly Panel B) was dropped
2026-10-05: both numbers it showed are already stated in Results 5.5 prose and in
Supplementary Table S12, so the panel was pure redundancy, not independent information.

Inputs: data/derived/intogen_rescoring/intogen_confidence_breakdown_detail.csv
Writes: figures/source/intogen_suppfig2_confidence_breakdown.pptx
"""
import pandas as pd
from pptx.enum.text import PP_ALIGN

from _style import DARK, GREY, LIGHT, M, OI, ROOT, line, new_slide, rect, save, text

SURF, MISS = OI["blue"], OI["vermillion"]


def main():
    der = ROOT / "data" / "derived" / "intogen_rescoring"
    det = pd.read_csv(der / "intogen_confidence_breakdown_detail.csv")

    prs, s = new_slide()
    text(s, M, M, 6.0, 0.25, "IntOGen's own confidence signals: surfaced vs never-surfaced driver pairs",
         pt=11, bold=True)
    x0, y0, w, h = M + 0.85, M + 0.65, 5.1, 2.3
    surf = det[det["surfaced"]]
    miss = det[~det["surfaced"]]
    surf_dist = (surf["n_methods"].value_counts(normalize=True).reindex(range(1, 8), fill_value=0) * 100)
    miss_dist = (miss["n_methods"].value_counts(normalize=True).reindex(range(1, 8), fill_value=0) * 100)

    for t in (0, 20, 40, 60, 80):
        yy = y0 + h - t / 80 * h
        line(s, x0, yy, x0 + w, yy, LIGHT, 0.5)
        text(s, x0 - 0.42, yy - 0.08, 0.4, 0.16, f"{t}%", pt=8, align=PP_ALIGN.RIGHT)
    line(s, x0, y0, x0, y0 + h, DARK, 0.75)
    line(s, x0, y0 + h, x0 + w, y0 + h, DARK, 0.75)

    n_cat = 7
    cw = w / n_cat
    bw = cw * 0.32
    for i, k in enumerate(range(1, 8)):
        cx = x0 + i * cw + cw / 2
        sv = surf_dist[k]
        mv = miss_dist[k]
        sh = sv / 80 * h
        mh = mv / 80 * h
        rect(s, cx - bw - 0.02, y0 + h - sh, bw, sh, SURF)
        rect(s, cx + 0.02, y0 + h - mh, bw, mh, MISS)
        text(s, cx - cw / 2, y0 + h + 0.04, cw, 0.16, str(k), pt=8, align=PP_ALIGN.CENTER)
    text(s, x0, y0 + h + 0.24, w, 0.16, "Number of IntOGen methods agreeing", pt=9, align=PP_ALIGN.CENTER)
    text(s, x0 - 0.42 - 0.1 - h / 2, y0 + h / 2 - 0.09, h, 0.18, "% of pairs", pt=9, align=PP_ALIGN.CENTER,
         rotate=270)

    lx, ly = x0 + w - 2.3, y0 + 0.05
    rect(s, lx, ly, 2.25, 0.5, "FFFFFF", line=GREY)
    for i, (lab, col) in enumerate(((f"Surfaced (n={len(surf)})", SURF), (f"Never surfaced (n={len(miss)})", MISS))):
        yy = ly + 0.13 + i * 0.22
        rect(s, lx + 0.08, yy, 0.18, 0.1, col)
        text(s, lx + 0.32, yy - 0.03, 1.8, 0.16, lab, pt=8, bold=True, colour=col)

    caption_y = y0 + h + 0.55
    text(s, M, caption_y, 6.0, 0.6,
         "Never-surfaced pairs skew toward IntOGen's single-method, lower-significance, "
         "lower-mutation-frequency tier (Mann-Whitney p<10⁻³⁵ for all three comparisons); "
         f"{100 * (miss['n_methods'] == 1).mean():.1f}% rest on only one of IntOGen's seven driver-detection methods, versus {100 * (surf['n_methods'] == 1).mean():.1f}% "
         "of surfaced pairs (Supplementary Table S12).", pt=9, italic=True, colour=DARK, wrap=True)

    save(prs, "intogen_suppfig2_confidence_breakdown.pptx")


if __name__ == "__main__":
    main()
