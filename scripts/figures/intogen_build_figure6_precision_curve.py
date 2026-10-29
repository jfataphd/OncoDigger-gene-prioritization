"""
intogen_build_figure6_precision_curve.py -- manuscript1.1 (IntOGen-primary draft) ONLY.
This is manuscript1.1's Figure 6, "IntOGen-driver enrichment across ranking depth", scored
against IntOGen's any-cancer driver set.

Inputs: data/derived/intogen_rescoring/precision_at_k_detail_intogen.csv,
        data/canonical/corpus_stats.csv,
        data/derived/intogen_rescoring/precision_at_k_bootstrap_ci_intogen.csv
Writes: figures/source/intogen_figure6_precision_curve.pptx
"""
import pandas as pd
from pptx.enum.dml import MSO_LINE_DASH_STYLE
from pptx.enum.text import PP_ALIGN

from _style import DARK, GREY, M, OI, ROOT, WHITE, axes, line, marker, new_slide, polygon, polyline, rect, save, text

import pandas as _pd
# distinct IntOGen driver genes (633) over the 19,250 protein-coding genes (Seal et al., 2026): the random baseline
_N_DRIVERS = _pd.read_csv(ROOT / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv",
                          sep=chr(9))["SYMBOL"].nunique()


TIERS = [(1, 10, "Ranks 1–10", "E8F1FA"), (10, 20, "Ranks 11–20", "E6F5EF"),
         (20, 30, "Ranks 21–30", "FDF3E1"), (30, 50, "Ranks 31–50", "F3EEF7")]


def main():
    der = ROOT / "data" / "derived" / "intogen_rescoring"
    ci = pd.read_csv(der / "precision_at_k_bootstrap_ci_intogen.csv").set_index("k")
    d = pd.read_csv(der / "precision_at_k_detail_intogen.csv")
    p10 = d[d["rank"] == 10].set_index("cancer")["precision_at_k"]
    p50 = d[d["rank"] == 50].set_index("cancer")["precision_at_k"]
    per_cancer = pd.DataFrame({"p10": p10, "p50": p50}).sort_values(["p10", "p50"], ascending=False)

    prs, s = new_slide()
    text(s, M, M, 6.0, 0.25, "Precision@k against IntOGen driver genes across 19 cancer types", pt=11, bold=True)
    x0, y0, w, h = M + 0.75, M + 0.75, 5.0, 4.3
    # tier shading first so everything else draws on top
    X = lambda k: x0 + (k - 1) / 49 * w
    for a, b, name, col in TIERS:
        rect(s, X(a), y0, X(b) - X(a), h, col)
        text(s, X(a), y0 - 0.2, X(b) - X(a), 0.16, name, pt=8, colour=GREY, align=PP_ALIGN.CENTER)
    P = axes(s, x0, y0, w, h, (1, 50), (0, 1), [1, 10, 20, 30, 40, 50], [0, 0.2, 0.4, 0.6, 0.8, 1.0],
             "OncoDigger gene rank (k)", "Cumulative IntOGen driver precision@k", fmt_y=lambda v: f"{v:.0%}", grid=False)

    ks = ci.index.tolist()
    polygon(s, [P(k, ci.loc[k, "ci_high_97.5"]) for k in ks] + [P(k, ci.loc[k, "ci_low_2.5"]) for k in ks[::-1]],
            OI["blue"], alpha=0.25)
    polyline(s, [P(k, ci.loc[k, "mean_p_at_k"]) for k in ks], "08306B", 2.25)
    ry = P(1, _N_DRIVERS / 19250)[1]
    line(s, x0, ry, x0 + w, ry, OI["vermillion"], 1.0, MSO_LINE_DASH_STYLE.DASH)

    # mean values at tier boundaries, labelled close to the curve, all at a consistent
    # offset below their own point; k=5 and k=10 are nudged apart horizontally (their
    # points are close together at this scale) rather than stacked at different heights.
    x_nudge = {5: -0.07, 10: 0.07}
    for k in (5, 10, 20, 30, 40, 50):
        px, py = P(k, ci.loc[k, "mean_p_at_k"])
        marker(s, px, py, 0.08, "08306B", line_colour=WHITE)
        label_y = py + 0.12
        text(s, px - 0.3 + x_nudge.get(k, 0.0), label_y, 0.6, 0.3,
             f"k={k}\n{ci.loc[k, 'mean_p_at_k']:.1%}", pt=7, align=PP_ALIGN.CENTER)
        line(s, px, py + 0.05, px, label_y - 0.01, GREY, 0.5, MSO_LINE_DASH_STYLE.ROUND_DOT)

    # per-cancer P@10/P@50 table, placed in the blank space below the CI band and above
    # the random-baseline line -- bounds computed from the actual lowest CI point and the
    # random line's own pixel position, not guessed, so it can never touch either.
    ci_low_max_y = max(P(k, ci.loc[k, "ci_low_2.5"])[1] for k in ks)
    table_top = ci_low_max_y + 0.35
    table_bottom = ry - 0.18
    n_left = 10
    left_rows, right_rows = per_cancer.iloc[:n_left], per_cancer.iloc[n_left:]
    n_rows = max(len(left_rows), len(right_rows)) + 1  # +1 for header
    row_h = (table_bottom - table_top) / n_rows
    block_w = w / 2 - 0.1
    name_w, num_w = block_w - 0.7, 0.35
    for bi, (rows, bx) in enumerate([(left_rows, x0 + 0.1), (right_rows, x0 + w / 2 + 0.1)]):
        yy = table_top
        text(s, bx, yy, name_w, row_h, "Cancer", pt=7, bold=True, colour=DARK)
        text(s, bx + name_w, yy, num_w, row_h, "P@10", pt=7, bold=True, colour=DARK, align=PP_ALIGN.CENTER)
        text(s, bx + name_w + num_w, yy, num_w, row_h, "P@50", pt=7, bold=True, colour=DARK, align=PP_ALIGN.CENTER)
        yy += row_h
        for cancer, r in rows.iterrows():
            text(s, bx, yy, name_w, row_h, cancer, pt=6.5, colour=DARK)
            text(s, bx + name_w, yy, num_w, row_h, f"{r['p10']:.0%}", pt=6.5, colour=DARK, align=PP_ALIGN.CENTER)
            text(s, bx + name_w + num_w, yy, num_w, row_h, f"{r['p50']:.0%}", pt=6.5, colour=DARK,
                 align=PP_ALIGN.CENTER)
            yy += row_h
        line(s, bx, table_top + row_h, bx + name_w + 2 * num_w, table_top + row_h, GREY, 0.5)

    # legend
    ly = y0 + h + 0.55
    line(s, x0, ly + 0.08, x0 + 0.3, ly + 0.08, "08306B", 2.25)
    text(s, x0 + 0.36, ly, 1.3, 0.16, "Mean (19 cancers)", pt=8)
    rect(s, x0 + 1.6, ly + 0.02, 0.25, 0.12, "BFD7EA")
    text(s, x0 + 1.9, ly, 1.3, 0.16, "95% bootstrap CI", pt=8)
    line(s, x0 + 3.1, ly + 0.08, x0 + 3.4, ly + 0.08, OI["vermillion"], 1.0, MSO_LINE_DASH_STYLE.DASH)
    text(s, x0 + 3.46, ly, 1.8, 0.16, f"Random ({_N_DRIVERS / 19250:.1%})", pt=8)
    save(prs, "intogen_figure6_precision_curve.pptx")


if __name__ == "__main__":
    main()
