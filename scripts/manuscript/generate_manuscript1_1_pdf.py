"""
generate_manuscript1_1_pdf.py -- manuscript1.1 (IntOGen-primary draft) ONLY. Builds a PDF
from OncoDigger_manuscript1.1.md, the comparison draft, using the same INSERT-marker-driven
rendering as generate_bioinformatics_submission_pdf.py (the real submission's PDF builder),
but pointed at manuscript1.1's own file and the figure PNGs (figures/png/intogen_*.png).

Output: manuscript/OncoDigger_manuscript1.1.pdf
Run:    python -X utf8 scripts/manuscript/generate_manuscript1_1_pdf.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, PageBreak, KeepTogether, HRFlowable,
    ListFlowable, ListItem,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import generate_pdf as _gp  # noqa: E402
from generate_pdf import (  # noqa: E402
    make_styles, md_inline_to_rml, parse_markdown, build_table_flowable, build_image_flowable,
    FIGS_DIR,
)

ROOT = Path(__file__).resolve().parents[2]
MD_PATH = ROOT / "manuscript" / "OncoDigger_manuscript1.1.md"
PDF_PATH = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "manuscript" / "OncoDigger_manuscript1.1.pdf"

PAGE_WIDTH, PAGE_HEIGHT = letter
MARGIN = 0.85 * inch
CONTENT_WIDTH = PAGE_WIDTH - 2 * MARGIN
_gp.CONTENT_WIDTH = CONTENT_WIDTH
_gp.TOP_MARGIN = MARGIN
_gp.BOTTOM_MARGIN = MARGIN

INSERT_RE = re.compile(r"<!--\s*INSERT (FIGURE \d+|TABLE \d+) HERE\s*-->")

# manuscript1.1's own figure-number -> IntOGen-rescored PNG mapping (parallel to
# generate_pdf.py's FIGURE_FILES / generate_bioinformatics_submission_pdf.py's
# SUPP_FIGURE_FILES_BIO, but pointing at the intogen_*.png files instead).
FIGURE_FILES_1_1 = {
    "Figure 1": ["intogen_figure1_workflow.png"],
    "Figure 2": ["intogen_figure2_gene_heatmap.png"],
    "Figure 3": ["intogen_figure3_validation.png"],
    "Figure 4": ["intogen_figure4_head_to_head.png"],
    "Figure 5": ["intogen_figure5_literature_weight.png"],
    "Figure 6": ["intogen_figure6_precision_curve.png"],
}
SUPP_FIGURE_FILES_1_1 = {
    "Supplementary Figure 1": ["intogen_suppfig1_pubmed_control.png"],
    "Supplementary Figure 2": ["intogen_suppfig2_confidence_breakdown.png"],
    "Supplementary Figure 3": ["intogen_suppfig3_recall_coverage.png"],
    "Supplementary Figure 4": ["intogen_suppfig4_retrieval_vocab_interaction.png"],
}


def parse_figure_legend_section(text: str) -> dict[str, str]:
    legends = {}
    blocks = re.split(r"\n(?=## )", text.strip())
    for block in blocks:
        block = block.strip()
        m = re.match(r"## (Figure \d+|Supplementary Figure S\d+)\.\s*([^\n]*)\n\n(.*)", block, re.DOTALL)
        if not m:
            continue
        label, title, body = m.group(1), m.group(2).strip(), m.group(3).strip()
        body = re.sub(r"\s*\n\s*---\s*$", "", body)  # a section divider after the last legend of a block is not part of it
        legends[label] = (title, body)
    return legends


def build_figure_flow(label: str, title_body: tuple[str, str], styles) -> list:
    files = FIGURE_FILES_1_1.get(label) or SUPP_FIGURE_FILES_1_1.get(
        label.replace("Supplementary Figure S", "Supplementary Figure ")
    ) or []
    flow = []
    for fname in files:
        path = FIGS_DIR / fname
        if path.exists():
            flow.extend(build_image_flowable(path, max_width=CONTENT_WIDTH))
            flow.append(Spacer(1, 4))
        else:
            print(f"  WARNING: missing image for {label}: {path}")
    if title_body:
        title, body = title_body
        flow.append(Paragraph(f"<b>{label}. {md_inline_to_rml(title)}.</b> {md_inline_to_rml(body)}",
                               styles["caption"]))
    return flow


def main() -> None:
    raw = MD_PATH.read_text(encoding="utf-8")
    raw = re.sub(r"<!--(?!\s*INSERT).*?-->\n?", "", raw, flags=re.DOTALL)

    # Captures ALL top-level "# Table N. ..." sections (there may be more than one, e.g.
    # Table 1 and Table 2), not just the first -- each runs until the next "# Table N"
    # heading or end of file, so adding a further table section never gets silently
    # swallowed into an earlier one's block.
    m_tables_all = re.search(r"\n# Table \d+\..*\Z", raw, re.DOTALL)
    table_blocks: dict[int, str] = {}
    if m_tables_all:
        tables_text = m_tables_all.group(0)
        raw = raw[: m_tables_all.start()]
        for m in re.finditer(r"^# Table (\d+)\.", tables_text, re.MULTILINE):
            n = int(m.group(1))
            start = m.start()
            next_m = re.search(r"^# Table \d+\.", tables_text[m.end():], re.MULTILINE)
            end = m.end() + next_m.start() if next_m else len(tables_text)
            table_blocks[n] = tables_text[start:end].strip()
    table1_block = table_blocks.get(1, "")

    m_legends = re.search(r"\n# Figure legends\n(.*?)(?=\n---\n\n# Table 1|\Z)", raw, re.DOTALL)
    legends = {}
    if m_legends:
        # Table 1 is already stripped out of `raw` above, so the "\n---\n\n# Table 1"
        # branch of the lookahead can never match here and this always falls through to
        # \Z (true end of the now-Table-1-less raw text) -- trimming a dangling trailing
        # "---" divider explicitly, since .strip() alone only removes whitespace, not "-".
        legend_text = re.sub(r"\n+---\s*\Z", "", m_legends.group(1))
        legends = parse_figure_legend_section(legend_text)
        raw = raw[: m_legends.start()]

    styles = make_styles()
    story = []

    title_match = re.search(r"^# (.+)$", raw, re.MULTILINE)
    main_title = title_match.group(1).strip()
    after_title = raw[title_match.end():]

    story.append(Paragraph(md_inline_to_rml(main_title), styles["title"]))

    byline_text = after_title.split("\n---\n", 1)[0]
    for block in parse_markdown(byline_text):
        if block[0] == "para":
            story.append(Paragraph(md_inline_to_rml(block[1]), styles["subtitle"]))
    story.append(HRFlowable(width="100%", thickness=0.5, color=colors.HexColor("#cbd5e0"),
                             spaceBefore=8, spaceAfter=8))

    remainder = after_title.split("\n---\n", 1)[1]

    lines = remainder.split("\n")
    buf: list[str] = []
    state = {"section": "", "ref_n": 0}

    def flush_buf():
        if not buf:
            return
        text = "\n".join(buf)
        buf.clear()
        for block in parse_markdown(text):
            kind = block[0]
            if kind == "h2":
                name = block[1]
                state["section"] = name.strip()
                if name.strip().upper() == "SUPPLEMENTARY MATERIAL":
                    story.append(PageBreak())
                    story.append(Paragraph("Supplementary Material", styles["title"]))
                else:
                    story.append(Paragraph(md_inline_to_rml(name), styles["h1"]))
            elif kind == "h3":
                story.append(Paragraph(md_inline_to_rml(block[1]), styles["h2"]))
            elif kind == "para":
                if state["section"] == "References":
                    story.append(Paragraph(md_inline_to_rml(block[1]), styles["ref"]))
                else:
                    story.append(Paragraph(md_inline_to_rml(block[1]), styles["body"]))
            elif kind == "list":
                items = [ListItem(Paragraph(md_inline_to_rml(t), styles["list_item"]), leftIndent=14)
                         for t in block[1]]
                story.append(ListFlowable(items, bulletType="bullet", bulletFontSize=8,
                                          bulletColor=colors.HexColor("#1565c0")))
                story.append(Spacer(1, 4))
            elif kind == "table":
                story.extend(build_table_flowable(block[1], styles))
            elif kind == "hrule":
                story.append(Spacer(1, 6))
                story.append(HRFlowable(width="100%", thickness=0.5, color=colors.HexColor("#cbd5e0")))
                story.append(Spacer(1, 6))

    for line in lines:
        m = INSERT_RE.search(line)
        if m:
            flush_buf()
            raw_label = m.group(1)
            label = raw_label.title()
            caption = legends.get(label)
            flow = build_figure_flow(label, caption, styles) if raw_label.startswith("FIGURE") else []
            m_table_label = re.match(r"TABLE (\d+)$", raw_label)
            if m_table_label and int(m_table_label.group(1)) in table_blocks:
                t_block = table_blocks[int(m_table_label.group(1))]
                t_title = re.search(r"^# (.+)$", t_block, re.MULTILINE).group(1)
                t_body = t_block[t_block.index("\n", 2):]
                flow = [Paragraph(md_inline_to_rml(t_title), styles["h2"])]
                for block in parse_markdown(t_body):
                    if block[0] == "table":
                        flow.extend(build_table_flowable(block[1], styles))
                    elif block[0] == "para":
                        flow.append(Paragraph(md_inline_to_rml(block[1]), styles["body"]))
            if flow:
                story.append(KeepTogether(flow))
                story.append(Spacer(1, 10))
            continue
        buf.append(line)
    flush_buf()

    supp_fig_labels = [l for l in legends if l.startswith("Supplementary Figure")]
    if supp_fig_labels:
        heading_pending = True
        for label in sorted(supp_fig_labels, key=lambda s: int(re.search(r"\d+", s).group())):
            flow = build_figure_flow(label, legends[label], styles)
            if flow:
                story.append(PageBreak())
                if heading_pending:
                    story.append(Paragraph("Supplementary Figures", styles["h1"]))
                    heading_pending = False
                story.extend(flow)  # no trailing Spacer: on a full page it would spill onto a blank page before the next PageBreak

    doc = SimpleDocTemplate(
        str(PDF_PATH), pagesize=letter,
        leftMargin=MARGIN, rightMargin=MARGIN,
        topMargin=MARGIN, bottomMargin=MARGIN,
        title=main_title, author="Jimmie Fata",
    )
    doc.build(story)
    print(f"Wrote: {PDF_PATH}")
    print(f"Size: {PDF_PATH.stat().st_size / 1024 / 1024:.2f} MB")


if __name__ == "__main__":
    main()
