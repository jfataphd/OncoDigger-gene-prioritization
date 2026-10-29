"""
export_pptx_to_png_highres.py
==============================
Export this repository's pptx figure sources to PNG at 300 DPI using PowerPoint's
Slide.Export with explicit pixel dimensions, not the default low-res SaveAs.

Computes pixel dimensions from the slide size (PowerPoint stores in points;
1 inch = 72 points; 300 DPI -> 300 px/in -> pixels = inches * 300).

Run:
  python -X utf8 scripts/manuscript/export_pptx_to_png_highres.py                  # all EXPORTS
  python -X utf8 scripts/manuscript/export_pptx_to_png_highres.py a.pptx b.pptx    # only these
"""
from __future__ import annotations

import sys
from pathlib import Path

import comtypes.client
from PIL import Image, ImageChops

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
PPTX_DIR = REPO_ROOT / "figures" / "source"
PNG_DIR = REPO_ROOT / "figures" / "png"

# Figures to (re-)export at high res. Add a filename here when its PPTX has been
# edited and the PNG needs regenerating.
EXPORTS = [
    "intogen_figure1_workflow.pptx",
    "intogen_figure2_gene_heatmap.pptx",
    "intogen_figure3_validation.pptx",
    "intogen_figure4_head_to_head.pptx",
    "intogen_figure5_literature_weight.pptx",
    "intogen_figure6_precision_curve.pptx",
    "intogen_suppfig1_pubmed_control.pptx",
    "intogen_suppfig2_confidence_breakdown.pptx",
    "intogen_suppfig3_recall_coverage.pptx",
    "intogen_suppfig4_retrieval_vocab_interaction.pptx",
]

TARGET_DPI = 300


def crop_bottom_whitespace(path: Path, padding: int = 60) -> None:
    """Trim trailing blank canvas from the bottom of an exported PNG (keeps full width
    and top margin untouched). Every slide is a fixed 6.7"x9.5" canvas regardless of how
    much content it actually holds, so a short single-panel figure still exports as a
    near-full-page-tall image; embedded at that height in the PDF, there is no room left
    on the page for its caption, which gets pushed onto an otherwise-blank following page.
    Cropping to actual content height lets the PDF scale the image to its real size and
    keep the caption on the same page."""
    img = Image.open(path).convert("RGB")
    bg = Image.new("RGB", img.size, (255, 255, 255))
    bbox = ImageChops.difference(img, bg).getbbox()
    if bbox is None:
        return
    content_bottom = min(img.height, bbox[3] + padding)
    if content_bottom < img.height - padding:
        img.crop((0, 0, img.width, content_bottom)).save(path)
        print(f"    Cropped bottom whitespace: {img.height}px -> {content_bottom}px")


def main():
    print("Starting PowerPoint COM...")
    ppt = comtypes.client.CreateObject("PowerPoint.Application")

    for fname in (sys.argv[1:] or EXPORTS):
        src = PPTX_DIR / fname
        if not src.exists():
            print(f"  MISSING: {src}")
            continue
        dst = PNG_DIR / (src.stem + ".png")

        print(f"  Opening: {src.name}")
        prs = ppt.Presentations.Open(str(src), WithWindow=False)

        slide_w_pt = prs.PageSetup.SlideWidth
        slide_h_pt = prs.PageSetup.SlideHeight
        slide_w_in = slide_w_pt / 72.0
        slide_h_in = slide_h_pt / 72.0
        px_w = int(round(slide_w_in * TARGET_DPI))
        px_h = int(round(slide_h_in * TARGET_DPI))
        print(f"    Slide: {slide_w_in:.2f}\" x {slide_h_in:.2f}\" -> {px_w} x {px_h} px @ {TARGET_DPI} DPI")

        slide = prs.Slides(1)  # 1-indexed in COM
        try:
            slide.Export(str(dst), "PNG", px_w, px_h)
            crop_bottom_whitespace(dst)
            print(f"    Saved: {dst.name} ({dst.stat().st_size // 1024} KB)")
        except Exception as e:
            print(f"    FAILED: {e}")

        prs.Close()

    ppt.Quit()
    print("Done.")


if __name__ == "__main__":
    main()
