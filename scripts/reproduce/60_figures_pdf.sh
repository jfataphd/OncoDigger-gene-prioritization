#!/bin/bash
# Stage 60: rebuild the figure sources (.pptx) from the regenerated data, then the manuscript PDF from the committed PNGs.
# This OVERWRITES figures/source/*.pptx: run it in a clone, not in a working copy whose PowerPoint files were edited by hand.
# PNG export (figures/png) needs Windows + Microsoft PowerPoint (scripts/manuscript/export_pptx_to_png_highres.py) and is NOT run here;
# without it the PDF uses the committed PNGs.
export REPRO_NO_CORPORA=1
. "$(dirname "$0")/env.sh" || exit 1
unset ONCODIGGER_ROOT BGFIX_ENABLE
for f in scripts/figures/intogen_build_*.py; do run f_$(basename "$f" .py) "$f"; done
run p01_pdf scripts/manuscript/generate_manuscript1_1_pdf.py
run p02_check_citations scripts/manuscript/check_citations.py
run p02b_check_references_exist scripts/manuscript/check_references_exist.py
run p03_check_text_numbers scripts/manuscript/check_text_numbers.py
run p04_check_provenance scripts/analysis/check_provenance.py
echo "STAGE 60 DONE"
