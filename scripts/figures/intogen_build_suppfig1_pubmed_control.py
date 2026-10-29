"""
intogen_build_suppfig1_pubmed_control.py -- manuscript1.1 (IntOGen-primary draft) ONLY.
This is manuscript1.1's Supplementary Figure S1, "cancer differentiation relative to the pooled
control". Panel A's cancer-specific segment is split by IntOGen any-cancer driver membership.
Panel B (Jaccard similarity) is pure gene-list overlap, independent of any reference standard.

Inputs: data/derived/intogen_rescoring/precision_at_k_detail_intogen.csv,
        data/derived/intogen_rescoring/corpus_specificity_clean_intogen.csv
Writes: figures/source/intogen_suppfig1_pubmed_control.pptx (6.7" x 9.5", BioSystems standard)

Run:  python -X utf8 scripts/figures/intogen_build_suppfig1_pubmed_control.py
"""
import itertools
from pathlib import Path

import pandas as pd
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Inches, Pt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures" / "source" / "intogen_suppfig1_pubmed_control.pptx"
_G = pd.read_csv(ROOT / "data" / "derived" / "pooled_control_mutation.csv")
GLOBAL_TOP10 = set(_G.loc[_G["rank"] <= 10, "gene"])   # pooled-corpus control (pooled_corpus_control.py)
SEG = [("shared", "Shared with pooled-corpus top-10", "BBBBBB"),
       ("specific_t1", "Cancer-specific, IntOGen driver", "0072B2"),
       ("specific_other", "Cancer-specific, not IntOGen driver", "E69F00")]
DARK, GREY = RGBColor(0x22, 0x22, 0x22), RGBColor(0x55, 0x55, 0x55)
LOW, HIGH = (0xFF, 0xFF, 0xFF), (0x08, 0x30, 0x6B)


def shade(v: float) -> RGBColor:
    t = max(0.0, min(1.0, v))
    return RGBColor(*[round(a + (b - a) * t) for a, b in zip(LOW, HIGH)])


def text(slide, x, y, w, h, s, pt=8, bold=False, colour=DARK, align=PP_ALIGN.LEFT,
         anchor=MSO_ANCHOR.MIDDLE, rot=0):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tb.rotation = rot
    tf = tb.text_frame
    tf.word_wrap = False
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    tf.vertical_anchor = anchor
    for i, line in enumerate(s.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        r = p.add_run()
        r.text = line
        r.font.name, r.font.size, r.font.bold, r.font.color.rgb = "Arial", Pt(pt), bold, colour
    return tb


def rect(slide, x, y, w, h, fill, line=None):
    r = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
    r.fill.solid()
    r.fill.fore_color.rgb = fill if isinstance(fill, RGBColor) else RGBColor.from_string(fill)
    if line is None:
        r.line.fill.background()
    else:
        r.line.color.rgb, r.line.width = line, Pt(0.5)
    r.shadow.inherit = False
    return r


def main() -> None:
    d = pd.read_csv(ROOT / "data" / "derived" / "intogen_rescoring" / "precision_at_k_detail_intogen.csv")
    top = d[d["rank"] <= 10]
    sets = top.groupby("cancer")["gene"].apply(set)
    intogen_any = top.groupby("cancer").apply(lambda g: set(g.loc[g["in_intogen_any"], "gene"]))
    comp = pd.DataFrame({
        "shared": sets.apply(lambda s: len(s & GLOBAL_TOP10)),
        "specific_t1": [len((sets[c] - GLOBAL_TOP10) & intogen_any[c]) for c in sets.index],
        "specific_other": [len((sets[c] - GLOBAL_TOP10) - intogen_any[c]) for c in sets.index],
    }, index=sets.index)
    check = pd.read_csv(ROOT / "data" / "derived" / "intogen_rescoring" / "corpus_specificity_clean_intogen.csv").set_index("cancer")
    assert (comp["shared"] == check.loc[comp.index, "overlap_with_global"]).all()
    assert (comp["specific_t1"] == check.loc[comp.index, "cancer_specific_intogen"]).all()
    comp = comp.sort_values(["shared", "specific_t1"], ascending=[True, False])

    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(6.7), Inches(9.5)
    s = prs.slides.add_slide(prs.slide_layouts[6])
    M = 0.35

    # ---------- (A) ----------
    text(s, M, M, 0.4, 0.25, "(A)", pt=12, bold=True)
    text(s, M + 0.42, M, 5.8, 0.25, "Composition of each OncoDigger top-10", pt=11, bold=True)
    lx, x0, x1, y0, bh, gap = M, 1.75, 5.55, M + 0.45, 0.13, 0.035
    for v in range(0, 11, 2):
        gx = x0 + (x1 - x0) * v / 10
        rect(s, gx - 0.003, y0 - 0.03, 0.006, len(comp) * (bh + gap), RGBColor(0xE5, 0xE5, 0xE5))
        text(s, gx - 0.2, y0 + len(comp) * (bh + gap) + 0.02, 0.4, 0.15, str(v), pt=8, colour=GREY, align=PP_ALIGN.CENTER)
    for i, (c, r) in enumerate(comp.iterrows()):
        y = y0 + i * (bh + gap)
        text(s, lx, y - 0.01, x0 - lx - 0.08, bh + 0.02, c, pt=8, align=PP_ALIGN.RIGHT)
        x = x0
        for key, _, col in SEG:
            w = (x1 - x0) * r[key] / 10
            if w > 0:
                rect(s, x, y, w, bh, col)
            x += w
    ybot = y0 + len(comp) * (bh + gap)
    text(s, x0, ybot + 0.2, x1 - x0, 0.18, "Genes in OncoDigger top-10", pt=10, align=PP_ALIGN.CENTER)
    for i, (_, lab, col) in enumerate(SEG):
        yy = ybot + 0.48 + i * 0.2
        rect(s, x0, yy + 0.02, 0.12, 0.12, col)
        text(s, x0 + 0.18, yy, 3.0, 0.16, lab, pt=8)

    # ---------- (B) ----------
    yB = ybot + 1.05
    text(s, M, yB, 0.4, 0.25, "(B)", pt=12, bold=True)
    text(s, M + 0.42, yB, 5.8, 0.25, "Similarity between cancer top-10 lists (Jaccard)", pt=11, bold=True)
    order = list(comp.index)
    short = [c.replace(" Cancer", "") for c in order]
    n = len(order)
    cell = 0.185
    hx0, hy0 = 1.35, yB + 0.4
    J = {}
    for a, b in itertools.product(order, order):
        J[a, b] = len(sets[a] & sets[b]) / len(sets[a] | sets[b])
    for i, a in enumerate(order):
        text(s, hx0 - 1.0, hy0 + i * cell, 0.93, cell, short[i], pt=7, align=PP_ALIGN.RIGHT)
        for j, b in enumerate(order):
            rect(s, hx0 + j * cell, hy0 + i * cell, cell, cell, shade(J[a, b]), line=RGBColor(0xFF, 0xFF, 0xFF))
    text(s, hx0, hy0 + n * cell + 0.05, n * cell, 0.15, "Columns in the same order as rows", pt=8,
         colour=GREY, align=PP_ALIGN.CENTER)
    off = [J[a, b] for a, b in itertools.combinations(order, 2)]
    mean = sum(off) / len(off)
    # colour scale + summary
    sx = hx0 + n * cell + 0.25
    steps = 40
    bar_h = n * cell * 0.45
    for k in range(steps):
        v = 1 - k / steps
        rect(s, sx, hy0 + bar_h * k / steps, 0.15, bar_h / steps + 0.002, shade(v))
    rect(s, sx, hy0, 0.15, bar_h, RGBColor(255, 255, 255), line=GREY).fill.background()
    text(s, sx + 0.2, hy0 - 0.06, 0.4, 0.14, "1.0", pt=8)
    text(s, sx + 0.2, hy0 + bar_h - 0.08, 0.4, 0.14, "0.0", pt=8)
    text(s, sx, hy0 + bar_h + 0.3, 1.1, 0.5,
         f"OncoDigger:\nmean Jaccard {mean:.3f}\nrange {min(off):.3f}–{max(off):.3f}", pt=8)
    text(s, sx, hy0 + bar_h + 1.0, 1.1, 0.35, "Pooled-corpus control:\n1.0 for every pair", pt=8)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    prs.save(OUT)
    print(f"Saved {OUT.relative_to(ROOT)}  mean Jaccard {mean:.3f} ({min(off):.3f}-{max(off):.3f}), "
          f"mean shared {comp['shared'].mean():.1f}")


if __name__ == "__main__":
    main()
