"""
intogen_build_suppfig4_retrieval_vocab_interaction.py -- manuscript1.1 (IntOGen-primary draft)
ONLY. New supplementary figure: per-cancer retrieval x vocabulary interaction term (Results 5.8,
Table 2), showing the interaction cancer by cancer (the number of cancers with a positive term is read from the
data, not typed). The axis extends below zero when any cancer's term is negative.

Inputs: data/derived/intogen_rescoring/interaction_all19_final.csv
Writes: figures/source/intogen_suppfig4_retrieval_vocab_interaction.pptx
"""
import pandas as pd
from pptx.enum.text import PP_ALIGN

from _style import DARK, LIGHT, M, OI, ROOT, line, new_slide, rect, save, text

BAR_COLOUR = OI["blue"]


def main():
    der = ROOT / "data" / "derived" / "intogen_rescoring"
    df = pd.read_csv(der / "interaction_all19_final.csv")
    df = df.sort_values("any_interaction", ascending=False).reset_index(drop=True)

    prs, s = new_slide()
    text(s, M, M, 6.0, 0.4,
         "Per-cancer retrieval x vocabulary interaction, any-cancer Precision@10",
         pt=11, bold=True, wrap=True)

    x0, y0, w, h = M + 0.55, M + 0.75, 5.75, 4.6
    ymax = max(8, int(df["any_interaction"].max()) + 1)
    ymin = min(0, int(df["any_interaction"].min()) - 1)
    span = ymax - ymin

    def ypos(v):
        return y0 + h - (v - ymin) / span * h

    for t in range(ymin, ymax + 1, 1):
        if t % 2 == 0:
            line(s, x0, ypos(t), x0 + w, ypos(t), LIGHT, 0.5)
            text(s, x0 - 0.45, ypos(t) - 0.08, 0.4, 0.16, f"+{t}" if t > 0 else f"{t}".replace("-", "−"), pt=8, align=PP_ALIGN.RIGHT)
    line(s, x0, y0, x0, y0 + h, DARK, 0.75)
    line(s, x0, ypos(0), x0 + w, ypos(0), DARK, 0.75)

    n = len(df)
    gw = w / n
    bw = gw * 0.62
    for i, row in df.iterrows():
        val = row["any_interaction"]
        xx = x0 + i * gw + (gw - bw) / 2
        top, bot = ypos(max(val, 0)), ypos(min(val, 0))
        if val != 0:
            rect(s, xx, top, bw, bot - top, BAR_COLOUR)
        lab = f"+{val}" if val > 0 else f"{val}".replace("-", "−")
        text(s, xx - 0.15, (top - 0.17) if val >= 0 else (bot + 0.01), bw + 0.3, 0.16, lab, pt=7, bold=True,
             colour=BAR_COLOUR, align=PP_ALIGN.CENTER)
        cx = x0 + i * gw + gw / 2
        label_w = 1.05
        cy = y0 + h + 0.05 + 0.55
        text(s, cx - label_w / 2, cy - gw / 2, label_w, gw, row["cancer"].replace(" Cancer", ""),
             pt=7, align=PP_ALIGN.RIGHT, rotate=270, wrap=False)

    text(s, x0 - 0.45 - 0.1 - h / 2, y0 + h / 2 - 0.09, h, 0.18,
         "Interaction (of 10 top-10 slots)", pt=10, align=PP_ALIGN.CENTER, rotate=270)

    save(prs, "intogen_suppfig4_retrieval_vocab_interaction.pptx")


if __name__ == "__main__":
    main()
