#!/bin/bash
# Stage 00: build the frozen corpus root (links to the 19 databases + vendored engine and lexicons) and check the corpus build.
. "$(dirname "$0")/env.sh" || exit 1
"$PYTHON" -X utf8 scripts/reproduce/assemble_roots.py --corpora "$ONCODIGGER_CORPORA" --work "$WORK" ${ONCODIGGER_RAW:+--raw "$ONCODIGGER_RAW"} || exit 1
run s00_check_corpus_build scripts/analysis/check_provenance.py
echo "STAGE 00 DONE"
