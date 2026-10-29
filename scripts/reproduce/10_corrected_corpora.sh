#!/bin/bash
# Stage 10: MET (Lung), REST (Thyroid) and FH (Kidney) corrected copies of three collections (about 4 GB, written to
# data/derived/intogen_rescoring/corrected_corpora, gitignored), then the corrected-collections root. The false-positive PMID sets come
# from the committed files in data/derived/intogen_rescoring/vocab_correction_pmids when the raw abstracts are absent, and are checked against
# them when ONCODIGGER_RAW is present.
. "$(dirname "$0")/env.sh" || exit 1
run s10_corrected_corpora scripts/analysis/vocab_corrected_rerun.py
"$PYTHON" -X utf8 scripts/reproduce/assemble_roots.py --corpora "$ONCODIGGER_CORPORA" --work "$WORK" --corr || exit 1
echo "STAGE 10 DONE"
