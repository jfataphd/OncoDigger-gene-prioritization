"""
generate_pdf.py — Build manuscript PDF from OncoDigger_reference.md
================================================================
Standard PDF output with figures embedded at first mention.

- Figures embedded as-is (no re-compression; ReportLab uses lossless Flate)
- Figures sized to page width preserving aspect ratio
- Tables rendered via ReportLab Table flowables
- Captions kept with figures via KeepTogether
- Tier 1 references rendered as numbered list

Output: manuscript/OncoDigger_reference.pdf (or the path given as the first argument)
Run:    python -X utf8 scripts/manuscript/generate_pdf.py [out.pdf]
"""
from __future__ import annotations

import re
from pathlib import Path

from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT, TA_JUSTIFY, TA_CENTER
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak,
    KeepTogether, Table, TableStyle, ListFlowable, ListItem,
)
from reportlab.platypus.flowables import HRFlowable

# === Paths ===
import sys

ROOT = Path(__file__).resolve().parent.parent.parent   # repo root
MD_PATH = ROOT / "manuscript" / "OncoDigger_reference.md"
# Optional first argument overrides the output path (e.g. for a test build)
PDF_PATH = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "manuscript" / "OncoDigger_reference.pdf"
FIGS_DIR = ROOT / "figures" / "png"

# === Figure -> file mapping (current numbering; see figures/README.md) ===
# Each main figure may have multiple panel PNGs; we embed them all stacked.
FIGURE_FILES = {
    "Figure 1": ["figure1_workflow.png"],
    "Figure 2": ["figure2_gene_heatmap.png"],
    "Figure 3": ["figure3_validation.png"],
    # NOTE: "Figure 4"/"Figure 5"/"Figure 6" are reader-facing labels, renumbered
    # so figures run in first-citation order; the underlying file names still
    # reflect their original content-based names, not their number. Swapped
    # 2026-09-30 (4/5) and again 2026-10-01 when the literature-weight figure
    # (Results 5.5) was inserted ahead of the ranking-depth figure.
    "Figure 4": ["figure5_head_to_head.png"],
    "Figure 5": ["figure6_literature_weight.png"],
    "Figure 6": ["figure4_precision_curve.png"],
}

# Candidates in order of preference (first existing file is used).
# Renumbered 2026-09-30 to run in this document's own citation order (was
# 3/4/6/7, inherited from a numbering pool shared across the 6-paper series;
# a standalone submission should not have unexplained gaps or out-of-order
# citations). The Bioinformatics doc uses its own separate numbering
# (SUPP_FIGURE_FILES_BIO in generate_bioinformatics_submission_pdf.py) since
# its condensed text cites the same figures in a different order.
SUPP_FIGURE_FILES = {
    "Supplementary Figure 1": ["suppfig3_precision_coverage.png"],
    "Supplementary Figure 2": ["suppfig6_pubmed_control.png"],
    "Supplementary Figure 3": ["suppfig7_chatgpt.png"],
    "Supplementary Figure 4": ["suppfig4_tier_precision.png"],
}

# === Page geometry ===
PAGE_WIDTH, PAGE_HEIGHT = letter
LEFT_MARGIN = 0.75 * inch
RIGHT_MARGIN = 0.75 * inch
TOP_MARGIN = 0.75 * inch
BOTTOM_MARGIN = 0.75 * inch
CONTENT_WIDTH = PAGE_WIDTH - LEFT_MARGIN - RIGHT_MARGIN  # ~7 inches


# === Styles ===
def make_styles():
    base = getSampleStyleSheet()
    styles = {}

    styles["title"] = ParagraphStyle(
        "title", parent=base["Title"],
        fontName="Helvetica-Bold", fontSize=18, leading=22,
        textColor=colors.HexColor("#0d1b2a"),
        alignment=TA_LEFT, spaceAfter=12,
    )
    styles["subtitle"] = ParagraphStyle(
        "subtitle", parent=base["Normal"],
        fontName="Helvetica", fontSize=11, leading=14,
        textColor=colors.HexColor("#4a5568"),
        spaceAfter=6,
    )
    styles["h1"] = ParagraphStyle(
        "h1", parent=base["Heading1"],
        fontName="Helvetica-Bold", fontSize=15, leading=18,
        textColor=colors.HexColor("#0d1b2a"),
        spaceBefore=18, spaceAfter=8, keepWithNext=True,
    )
    styles["h2"] = ParagraphStyle(
        "h2", parent=base["Heading2"],
        fontName="Helvetica-Bold", fontSize=12, leading=15,
        textColor=colors.HexColor("#1a365d"),
        spaceBefore=14, spaceAfter=6, keepWithNext=True,
    )
    styles["h3"] = ParagraphStyle(
        "h3", parent=base["Heading3"],
        fontName="Helvetica-Bold", fontSize=10.5, leading=13,
        textColor=colors.HexColor("#2d3748"),
        spaceBefore=8, spaceAfter=4, keepWithNext=True,
    )
    styles["body"] = ParagraphStyle(
        "body", parent=base["BodyText"],
        fontName="Times-Roman", fontSize=10.5, leading=14,
        alignment=TA_JUSTIFY,
        spaceAfter=8,
    )
    styles["abstract"] = ParagraphStyle(
        "abstract", parent=base["BodyText"],
        fontName="Times-Roman", fontSize=10, leading=13,
        alignment=TA_JUSTIFY, textColor=colors.HexColor("#2d3748"),
        spaceAfter=8, leftIndent=12, rightIndent=12,
    )
    styles["caption"] = ParagraphStyle(
        "caption", parent=base["BodyText"],
        fontName="Helvetica", fontSize=8.5, leading=11.5,
        alignment=TA_JUSTIFY, textColor=colors.HexColor("#4a5568"),
        spaceBefore=4, spaceAfter=14, leftIndent=4, rightIndent=4,
    )
    styles["fig_label"] = ParagraphStyle(
        "fig_label", parent=base["BodyText"],
        fontName="Helvetica-Bold", fontSize=9, leading=12,
        textColor=colors.HexColor("#0d1b2a"), spaceBefore=4, spaceAfter=2,
    )
    styles["ref"] = ParagraphStyle(
        "ref", parent=base["BodyText"],
        fontName="Helvetica", fontSize=9, leading=12,
        alignment=TA_LEFT, leftIndent=18, firstLineIndent=-18,
        spaceAfter=3,
    )
    styles["list_item"] = ParagraphStyle(
        "list_item", parent=base["BodyText"],
        fontName="Helvetica", fontSize=10, leading=13,
        leftIndent=14, bulletIndent=2,
        spaceAfter=3,
    )
    styles["highlight"] = ParagraphStyle(
        "highlight", parent=base["BodyText"],
        fontName="Helvetica", fontSize=10, leading=13,
        textColor=colors.HexColor("#1a365d"),
        leftIndent=14, bulletIndent=2,
        spaceAfter=3,
    )
    return styles


# === Markdown inline formatting -> ReportLab markup ===
def md_inline_to_rml(text: str) -> str:
    """Convert markdown inline to ReportLab Paragraph XML-like markup."""
    # Escape XML-special chars first (but leave existing entities)
    text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    # Bold + italic
    text = re.sub(r"\*\*\*(.+?)\*\*\*", r"<b><i>\1</i></b>", text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<!\w)\*([^*\n]+?)\*(?!\w)", r"<i>\1</i>", text)
    # Code spans -> monospace
    text = re.sub(r"`([^`]+)`", r"<font face='Courier' size='9'>\1</font>", text)
    # Superscript. Exponents are written as a single trailing caret (e.g.
    # "10^-8"), not a matched pair, so match just the caret plus the numeric
    # exponent -- a paired-caret pattern here would instead span from one
    # exponent's caret to the next one's, swallowing everything between them
    # into a single stray superscript run. No subscript markup is used
    # anywhere in this manuscript's prose, and a paired-tilde "~x~" pattern
    # (the same risk class) caused exactly that bug with "~" used as a plain
    # "approximately" symbol, so there is no subscript rule here.
    text = re.sub(r"\^(-?\d+)", r"<super>\1</super>", text)
    return text


# === Parse markdown into ordered blocks ===
def parse_markdown(text: str):
    """Yield blocks of type: ('section_h1', name), ('section_h2', name), ('section_h3', name),
    ('para', text), ('list', items), ('table', rows), ('hrule',)"""
    lines = text.split("\n")
    i = 0
    in_ul = False
    list_buf: list[str] = []
    para_buf: list[str] = []
    table_buf: list[list[str]] = []
    in_table = False

    def flush_para():
        nonlocal para_buf
        if para_buf:
            joined = " ".join(l.strip() for l in para_buf if l.strip())
            if joined:
                yield_buf.append(("para", joined))
            para_buf = []

    def flush_list():
        nonlocal in_ul, list_buf
        if list_buf:
            yield_buf.append(("list", list_buf))
            list_buf = []
        in_ul = False

    def flush_table():
        nonlocal in_table, table_buf
        if table_buf:
            yield_buf.append(("table", table_buf))
            table_buf = []
        in_table = False

    yield_buf: list = []

    for line in lines:
        s = line.rstrip()
        stripped = s.strip()

        # Section headers
        m1 = re.match(r"^# (.+)$", stripped)
        m2 = re.match(r"^## (.+)$", stripped)
        m3 = re.match(r"^### (.+)$", stripped)
        if m1:
            flush_para()
            flush_list()
            flush_table()
            yield_buf.append(("title", m1.group(1)))
            continue
        if m2:
            flush_para()
            flush_list()
            flush_table()
            yield_buf.append(("h2", m2.group(1)))
            continue
        if m3:
            flush_para()
            flush_list()
            flush_table()
            yield_buf.append(("h3", m3.group(1)))
            continue

        # Horizontal rule
        if re.match(r"^-{3,}$", stripped):
            flush_para()
            flush_list()
            flush_table()
            yield_buf.append(("hrule",))
            continue

        # Table row
        if stripped.startswith("|") and stripped.endswith("|") and stripped.count("|") >= 2:
            flush_para()
            flush_list()
            in_table = True
            cells = [c.strip() for c in stripped.strip("|").split("|")]
            table_buf.append(cells)
            continue
        if in_table and not (stripped.startswith("|") and stripped.endswith("|")):
            flush_table()

        # List items
        if stripped.startswith("- ") or stripped.startswith("* "):
            flush_para()
            in_ul = True
            list_buf.append(stripped[2:].strip())
            continue
        if in_ul and not stripped:
            flush_list()
            continue
        if in_ul and stripped and not stripped.startswith(("- ", "* ")):
            flush_list()

        # Blank line -> end paragraph
        if not stripped:
            flush_para()
            continue

        # Regular paragraph line
        para_buf.append(stripped)

    flush_para()
    flush_list()
    flush_table()
    return yield_buf


# === Build flowables ===
def build_image_flowable(path: Path, caption_label: str = "", max_width: float = None):
    """Return list of flowables: Image + caption label (if provided)."""
    if not path.exists():
        return []
    if max_width is None:
        max_width = CONTENT_WIDTH
    from reportlab.lib.utils import ImageReader
    ir = ImageReader(str(path))
    iw, ih = ir.getSize()
    # Scale to fit max_width while preserving aspect ratio
    scale = min(max_width / iw, (PAGE_HEIGHT - TOP_MARGIN - BOTTOM_MARGIN - 2.1*inch) / ih)
    if scale > 1.0:
        scale = 1.0
    w = iw * scale
    h = ih * scale
    img = Image(str(path), width=w, height=h)
    return [img]


def build_figure_block(fig_label: str, caption_html: str, styles):
    """Return KeepTogether flowable for a figure: images + caption + label."""
    flow = []
    img_paths = FIGURE_FILES.get(fig_label, [])
    for fname in img_paths:
        path = FIGS_DIR / fname
        if path.exists():
            flow.extend(build_image_flowable(path))
            flow.append(Spacer(1, 4))
    if not flow:
        return []
    # Caption
    flow.append(Paragraph(f"<b>{fig_label}.</b> {caption_html}", styles["caption"]))
    return [KeepTogether(flow), Spacer(1, 6)]


def build_supp_figure_block(fig_label: str, caption_html: str, styles):
    flow = []
    for fname in SUPP_FIGURE_FILES.get(fig_label, []):
        path = FIGS_DIR / fname
        if path.exists():
            flow.extend(build_image_flowable(path))
            flow.append(Spacer(1, 4))
            break
    else:
        print(f"  WARNING: no image found for {fig_label}")
    flow.append(Paragraph(f"<b>{fig_label}.</b> {caption_html}", styles["caption"]))
    return [KeepTogether(flow), Spacer(1, 6)]


def build_table_flowable(rows: list[list[str]], styles):
    """Render a markdown table as a ReportLab Table.

    Convention: first row = header, second row = separator (skipped),
    remaining = body rows.
    """
    if len(rows) < 2:
        return []
    header = rows[0]
    body = rows[2:] if len(rows) > 1 else []

    # Convert each cell via md_inline_to_rml; wrap in Paragraphs for wrapping
    def cell(text: str, header: bool = False) -> Paragraph:
        style = ParagraphStyle(
            "tcell", fontName="Helvetica-Bold" if header else "Helvetica",
            fontSize=8.5, leading=11,
            textColor=colors.HexColor("#0d1b2a") if header else colors.HexColor("#2d3748"),
        )
        return Paragraph(md_inline_to_rml(text), style)

    data = [[cell(c, True) for c in header]]
    for row in body:
        # Pad row to header length
        padded = list(row) + [""] * max(0, len(header) - len(row))
        data.append([cell(c, False) for c in padded[: len(header)]])

    n_cols = len(header)
    col_w = CONTENT_WIDTH / n_cols
    tbl = Table(data, colWidths=[col_w] * n_cols, repeatRows=1)
    tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#edf2f7")),
        ("LINEBELOW", (0, 0), (-1, 0), 1.2, colors.HexColor("#a0aec0")),
        ("LINEBELOW", (0, 1), (-1, -1), 0.3, colors.HexColor("#e2e8f0")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1),
         [colors.white, colors.HexColor("#fafbfc")]),
    ]))
    return [tbl, Spacer(1, 8)]


# === Parse figure legends section into {label: caption} ===
def parse_figure_legends(text: str) -> dict[str, str]:
    legends = {}
    blocks = re.split(r"\n(?=\*\*Figure\s)", text.strip())
    for block in blocks:
        block = block.strip()
        m = re.match(r"\*\*(Figure \d+)\.([^*]*)\*\*\s*(.*)", block, re.DOTALL)
        if not m:
            continue
        label = m.group(1).strip()
        title = m.group(2).strip()
        body = m.group(3).strip()
        full = f"<b>{md_inline_to_rml(title)}</b> {md_inline_to_rml(body)}"
        legends[label] = full
    return legends


def parse_supp_figure_legends(text: str) -> tuple[dict[str, str], list[tuple[str, str]]]:
    """Returns (supp_fig_captions, supp_table_blocks).
    supp_table_blocks is a list of (label, body_markdown_text) preserving order.
    """
    captions = {}
    table_blocks = []
    blocks = re.split(r"\n(?=\*\*Supplementary (?:Figure|Table)\s)", text.strip())
    for block in blocks:
        block = block.strip()
        m_fig = re.match(r"\*\*(Supplementary Figure \d+)\.\*\*\s*(.*)", block, re.DOTALL)
        m_tbl = re.match(r"\*\*(Supplementary Table \d+)\.\*\*\s*(.*)", block, re.DOTALL)
        if m_fig:
            captions[m_fig.group(1)] = md_inline_to_rml(m_fig.group(2).strip())
        elif m_tbl:
            table_blocks.append((m_tbl.group(1), m_tbl.group(2).strip()))
    return captions, table_blocks


# === Main build ===
def main():
    md_text = MD_PATH.read_text(encoding="utf-8")
    styles = make_styles()

    # Split into top-level sections (## Headers)
    section_pattern = re.compile(r"^## (.+)$", re.MULTILINE)
    sections = {}
    matches = list(section_pattern.finditer(md_text))
    if matches:
        for i, m in enumerate(matches):
            start = m.end()
            end = matches[i + 1].start() if i + 1 < len(matches) else len(md_text)
            name = re.sub(r"^\d+\.\s+", "", m.group(1).strip())
            # Drop "---" section dividers: they are Markdown separators, not content
            sections[name] = re.sub(r"(?m)^---\s*$", "", md_text[start:end])

    # Title (the # title line at the top)
    title_match = re.search(r"^# (.+)$", md_text, re.MULTILINE)
    main_title = title_match.group(1).strip() if title_match else "Manuscript"

    # Figure legends parsed once for use in inline figure insertion
    fig_legends = parse_figure_legends(sections.get("Figure Legends", ""))

    # Build flowables
    story = []
    referenced_figs: set[str] = set()

    def insert_figures_for_para(text: str):
        """If the paragraph mentions a Figure N not yet inserted, add it."""
        out = []
        # Find all Figure N references (also "Figure NA" → main figure)
        refs = re.findall(r"(?<!Supplementary )\bFigures? (\d+)[A-Z]?\b", text)
        for ref_num in refs:
            fig_label = f"Figure {ref_num}"
            if fig_label in referenced_figs:
                continue
            if fig_label not in FIGURE_FILES:
                continue
            caption = fig_legends.get(fig_label, "")
            block = build_figure_block(fig_label, caption, styles)
            out.extend(block)
            referenced_figs.add(fig_label)
        return out

    def render_section(name: str, text: str, section_title: str | None = None):
        if section_title:
            story.append(Paragraph(section_title, styles["h1"]))
        blocks = parse_markdown(text)
        for kind, *payload in (b if isinstance(b, tuple) else (b,) for b in blocks):
            if kind == "h2":
                story.append(Paragraph(md_inline_to_rml(payload[0]), styles["h2"]))
            elif kind == "h3":
                story.append(Paragraph(md_inline_to_rml(payload[0]), styles["h3"]))
            elif kind == "para":
                story.append(Paragraph(md_inline_to_rml(payload[0]), styles["body"]))
                # Insert any newly-mentioned figures after this paragraph
                story.extend(insert_figures_for_para(payload[0]))
            elif kind == "list":
                items = [
                    ListItem(Paragraph(md_inline_to_rml(it), styles["list_item"]), leftIndent=14)
                    for it in payload[0]
                ]
                story.append(ListFlowable(items, bulletType="bullet", bulletFontSize=8,
                                          bulletColor=colors.HexColor("#1565c0")))
                story.append(Spacer(1, 4))
            elif kind == "table":
                story.extend(build_table_flowable(payload[0], styles))
            elif kind == "hrule":
                story.append(Spacer(1, 6))
                story.append(HRFlowable(width="100%", thickness=0.5,
                                         color=colors.HexColor("#cbd5e0")))
                story.append(Spacer(1, 6))

    # === Title block ===
    story.append(Paragraph(main_title, styles["title"]))
    story.append(Paragraph("<b>Jimmie Fata</b>", styles["subtitle"]))
    story.append(Paragraph(
        "Department of Biology, College of Staten Island, "
        "City University of New York, Staten Island, New York, USA",
        styles["subtitle"]))
    story.append(Paragraph("Corresponding author: jimmie.fata@csi.cuny.edu",
                            styles["subtitle"]))
    story.append(HRFlowable(width="100%", thickness=0.5,
                             color=colors.HexColor("#cbd5e0"),
                             spaceBefore=8, spaceAfter=8))

    # === Abstract ===
    abstract_text = sections.get("Abstract", "").strip()
    if abstract_text:
        story.append(Paragraph("Abstract", styles["h1"]))
        # Render as one block paragraph + Keywords if separated
        # The abstract is mostly one paragraph followed by **Keywords:** ...
        # Split on **Keywords:** if present
        kw_match = re.search(r"\*\*Keywords:\*\*\s*(.*)", abstract_text, re.DOTALL)
        body_only = abstract_text
        kw_only = None
        if kw_match:
            body_only = abstract_text[:kw_match.start()].strip()
            kw_only = kw_match.group(1).strip()
        for para in re.split(r"\n\s*\n", body_only):
            para = para.strip()
            if para:
                story.append(Paragraph(md_inline_to_rml(para), styles["abstract"]))
        if kw_only:
            story.append(Paragraph(f"<b>Keywords:</b> {md_inline_to_rml(kw_only)}",
                                    styles["abstract"]))

    # === Highlights ===
    hl_text = sections.get("Highlights", "").strip()
    if hl_text:
        story.append(Paragraph("Highlights", styles["h1"]))
        # Each bullet is a list item
        items = []
        for line in hl_text.split("\n"):
            line = line.strip()
            if line.startswith("- ") or line.startswith("* "):
                items.append(ListItem(
                    Paragraph(md_inline_to_rml(line[2:].strip()), styles["highlight"]),
                    leftIndent=14))
            elif line.startswith("*(") and line.endswith(")*"):
                # Italic note line
                story.append(Paragraph(f"<i>{md_inline_to_rml(line[1:-1])}</i>", styles["abstract"]))
        if items:
            story.append(ListFlowable(items, bulletType="bullet", bulletFontSize=9,
                                       bulletColor=colors.HexColor("#1565c0")))

    # === Main sections (Intro, Methods, Results, Discussion) ===
    for section_name in ["Introduction", "Methods", "Results", "Discussion"]:
        sec_text = sections.get(section_name, "")
        if sec_text:
            render_section(section_name, sec_text, section_title=section_name)

    # === Insert any figures not yet referenced (rare, e.g., if not mentioned by number) ===
    # Skipped intentionally; if a figure isn't referenced in text, it won't be inserted inline.

    # === Figure Legends ===
    legend_text = sections.get("Figure Legends", "")
    if legend_text:
        story.append(PageBreak())
        story.append(Paragraph("Figure Legends", styles["h1"]))
        # Each legend block: **Figure N. Title** body...
        blocks = re.split(r"\n(?=\*\*Figure\s)", legend_text.strip())
        for block in blocks:
            block = block.strip()
            if not block:
                continue
            m = re.match(r"\*\*(Figure \d+)\.([^*]*)\*\*\s*(.*)", block, re.DOTALL)
            if not m:
                continue
            label = m.group(1).strip()
            title = m.group(2).strip()
            body = m.group(3).strip()
            story.append(Paragraph(
                f"<b>{label}. {md_inline_to_rml(title)}</b> {md_inline_to_rml(body)}",
                styles["body"]))

    # === Prose sections ===
    for section_name in [
        "CRediT Author Contribution Statement",
        "Declaration of Competing Interests",
        "Funding",
        "Declaration of Generative AI and AI-Assisted Technologies in the Manuscript Preparation Process",
        "Data and Code Availability Statement",
        "Acknowledgements",
    ]:
        sec_text = sections.get(section_name, "").strip()
        if sec_text:
            story.append(Paragraph(section_name, styles["h2"]))
            for para in re.split(r"\n\s*\n", sec_text):
                para = para.strip()
                if para:
                    story.append(Paragraph(md_inline_to_rml(para), styles["body"]))

    # === References ===
    ref_text = sections.get("References", "").strip()
    if ref_text:
        story.append(PageBreak())
        story.append(Paragraph("References", styles["h1"]))
        ref_paras = [p.strip() for p in re.split(r"\n\s*\n", ref_text) if p.strip()]
        for i, ref in enumerate(ref_paras, 1):
            joined = " ".join(l.strip() for l in ref.splitlines() if l.strip())
            story.append(Paragraph(f"{i}. {md_inline_to_rml(joined)}", styles["ref"]))

    # === Supplementary Figures + Tables ===
    supp_text = sections.get("Supplementary Figures", "").strip()
    if supp_text:
        story.append(PageBreak())
        story.append(Paragraph("Supplementary Materials", styles["h1"]))

        # Walk block-by-block in original order: figures and tables interleaved
        blocks = re.split(r"\n(?=\*\*Supplementary (?:Figure|Table)\s)", supp_text)
        for block in blocks:
            block = block.strip()
            if not block:
                continue
            m_fig = re.match(r"\*\*(Supplementary Figure \d+)\.\*\*\s*(.*)", block, re.DOTALL)
            m_tbl = re.match(r"\*\*(Supplementary Table \d+)\.\*\*\s*(.*)", block, re.DOTALL)
            if m_fig:
                label = m_fig.group(1)
                caption = md_inline_to_rml(m_fig.group(2).strip())
                story.extend(build_supp_figure_block(label, caption, styles))
            elif m_tbl:
                label = m_tbl.group(1)
                body = m_tbl.group(2).strip()
                # The body may contain caption text + a markdown table
                # Find the first table row line, render caption above, table below
                lines = body.split("\n")
                cap_lines = []
                tbl_start = None
                for j, ln in enumerate(lines):
                    if ln.strip().startswith("|"):
                        tbl_start = j
                        break
                    cap_lines.append(ln)
                # Everything after the table is treated as note paragraph
                cap_text = " ".join(l.strip() for l in cap_lines if l.strip())
                story.append(Paragraph(
                    f"<b>{label}.</b> {md_inline_to_rml(cap_text)}",
                    styles["caption"]))
                if tbl_start is not None:
                    # Collect contiguous table lines
                    tbl_lines = []
                    j = tbl_start
                    while j < len(lines) and lines[j].strip().startswith("|"):
                        cells = [c.strip() for c in lines[j].strip().strip("|").split("|")]
                        tbl_lines.append(cells)
                        j += 1
                    if tbl_lines:
                        story.extend(build_table_flowable(tbl_lines, styles))
                    # After-table notes
                    remaining = " ".join(
                        l.strip() for l in lines[j:] if l.strip()
                    )
                    if remaining:
                        story.append(Paragraph(md_inline_to_rml(remaining), styles["caption"]))

    # === Render PDF ===
    doc = SimpleDocTemplate(
        str(PDF_PATH), pagesize=letter,
        leftMargin=LEFT_MARGIN, rightMargin=RIGHT_MARGIN,
        topMargin=TOP_MARGIN, bottomMargin=BOTTOM_MARGIN,
        title=main_title, author="Jimmie Fata",
    )
    doc.build(story)
    print(f"Wrote: {PDF_PATH}")
    print(f"Size: {PDF_PATH.stat().st_size / 1024 / 1024:.2f} MB")


if __name__ == "__main__":
    main()
