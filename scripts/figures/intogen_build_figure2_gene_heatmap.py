"""
intogen_build_figure2_gene_heatmap.py -- manuscript1.1 (IntOGen-primary draft) ONLY.
The gene/rank grid of the top-50 rankings (the ranking itself is reference-agnostic); the strip
marks IntOGen any-cancer driver membership.

Input:  data/derived/intogen_rescoring/precision_at_k_detail_intogen.csv, data/canonical/corpus_stats.csv
Writes: figures/source/intogen_figure2_gene_heatmap.pptx
"""
import pandas as pd
from pptx.enum.dml import MSO_LINE_DASH_STYLE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN

from _style import (DARK, GREY, M, OI, ROOT, SYSTEM_COLOURS, WHITE, blues, line, new_slide, rect, rgb, save, text)

ORDER = ["Hematologic", "Gynecologic/Breast", "Gastrointestinal", "Genitourinary", "Thoracic", "Skin", "Endocrine", "CNS"]


def shade(r: float):
    return WHITE if pd.isna(r) else blues(0.15 + 0.85 * (20 - r) / 19)


def main():
    d = pd.read_csv(ROOT / "data" / "derived" / "intogen_rescoring" / "precision_at_k_detail_intogen.csv")
    system = pd.read_csv(ROOT / "data" / "canonical" / "corpus_stats.csv").set_index("cancer")["system"]
    top = d[d["rank"] <= 20]
    rk = top.pivot(index="cancer", columns="gene", values="rank")
    intogen_any = set(d.loc[d["in_intogen_any"], "gene"])
    n = rk.notna().sum()
    best = rk.min()
    genes = n[n >= 2].index
    est = sorted([g for g in genes if best[g] <= 10], key=lambda g: (-n[g], best[g], g))
    emg = sorted([g for g in genes if best[g] > 10], key=lambda g: (-n[g], best[g], g))
    cols = est + emg
    rows = sorted(rk.index, key=lambda c: (ORDER.index(system[c]), c))

    prs, s = new_slide()
    text(s, M, M, 6.0, 0.25, "Pan-cancer top-20 mutation-associated gene signatures", pt=11, bold=True)
    lab_w = 1.3
    x0 = M + lab_w + 0.05
    grid_w = 6.7 - M - 0.55 - x0
    cw = grid_w / len(cols)
    ch = 0.28
    y_tier = M + 0.75
    y0 = y_tier + 0.22

    # tier blocks
    xe = x0 + len(est) * cw
    rect(s, x0, M + 0.38, len(est) * cw, 0.3, "DCEBF7", line=OI["sky"])
    text(s, x0, M + 0.38, len(est) * cw, 0.3, "Genes in the top 10 of ≥1 cancer", pt=8, bold=True,
         colour=OI["blue"], align=PP_ALIGN.CENTER)
    if emg:
        rect(s, xe, M + 0.38, len(emg) * cw, 0.3, "D5F2E8", line=OI["green"])
        text(s, xe, M + 0.38, len(emg) * cw, 0.3, "11–20\nonly", pt=7, bold=True, colour="00664A",
             align=PP_ALIGN.CENTER)
    # IntOGen driver strip
    for j, g in enumerate(cols):
        rect(s, x0 + j * cw, y_tier, cw, 0.16, "333333" if g in intogen_any else "DDDDDD", line=WHITE, lw=0.25)
    text(s, M, y_tier - 0.02, lab_w, 0.2, "IntOGen driver", pt=8, align=PP_ALIGN.RIGHT)

    # cells
    for i, c in enumerate(rows):
        y = y0 + i * ch
        rect(s, M, y + 0.01, lab_w, ch - 0.02, SYSTEM_COLOURS[system[c]])
        text(s, M, y, lab_w - 0.06, ch, c, pt=8, bold=True, colour=WHITE, align=PP_ALIGN.RIGHT)
        for j, g in enumerate(cols):
            rect(s, x0 + j * cw, y, cw, ch, shade(rk.loc[c, g] if g in rk.columns else None), line="BBBBBB", lw=0.25)
    ybot = y0 + len(rows) * ch
    rect(s, x0, y0, grid_w, ybot - y0, None, line=DARK, lw=0.75)
    if emg:
        line(s, xe, M + 0.38, xe, ybot, OI["vermillion"], 1.0, MSO_LINE_DASH_STYLE.DASH)

    # gene labels (rotated)
    half = len(rows) / 2
    # staggered on two levels so adjacent 7-pt labels do not collide
    for j, g in enumerate(cols):
        off = 0.0 if j % 2 == 0 else 0.62
        cx, cy = x0 + j * cw + cw / 2, ybot + 0.04 + off + 0.3
        if off:
            line(s, cx, ybot, cx, ybot + off, "BBBBBB", 0.25)
        text(s, cx - 0.3, cy - cw / 2 - 0.01, 0.6, cw + 0.02, g, pt=7, bold=n[g] >= half,
             align=PP_ALIGN.RIGHT, rotate=270, anchor=MSO_ANCHOR.MIDDLE)

    # rank scale
    lx = 6.7 - M - 0.4
    text(s, lx - 0.05, y0, 0.5, 0.16, "Rank", pt=8, bold=True)
    for k, r in enumerate([1, 5, 10, 15, 20, None]):
        yy = y0 + 0.25 + k * 0.3
        rect(s, lx, yy, 0.16, 0.22, shade(r), line=GREY, lw=0.5)
        text(s, lx + 0.2, yy, 0.3, 0.22, str(r) if r else "—", pt=8)
    text(s, lx - 0.05, y0 + 0.25 + 6 * 0.3, 0.5, 0.3, "— not in\ntop-20", pt=7, colour=GREY)

    # organ-system legend
    ly = ybot + 1.45
    text(s, M, ly, 1.3, 0.16, "Organ system", pt=8, bold=True)
    for k, sys_ in enumerate(ORDER):
        x = M + (k % 4) * 1.55
        yy = ly + 0.22 + (k // 4) * 0.2
        rect(s, x, yy + 0.02, 0.12, 0.12, SYSTEM_COLOURS[sys_])
        text(s, x + 0.17, yy, 1.35, 0.16, sys_, pt=8)
    text(s, M, ly + 0.7, 6.0, 0.48,
         "Column blocks group genes; cell shading gives each cancer's own rank. Right block: genes reaching\n"
         "ranks 11–20 only. IntOGen driver strip: dark = gene listed by IntOGen as a driver (any cancer type),\n"
         "light grey = not listed. Bold gene symbols: in the top-20 of at least half of the 19 cancers.", pt=8,
         colour=DARK, wrap=True)
    save(prs, "intogen_figure2_gene_heatmap.pptx")
    print(f"{len(est)} established + {len(emg)} Tier 2 columns; bold: {[g for g in cols if n[g] >= half]}")


if __name__ == "__main__":
    main()
