#!/bin/bash
# Run every stage in order. Needs ONCODIGGER_CORPORA (see env.sh). Stages 00-40 need the corpus databases; 50-60 do not.
set -e
d="$(dirname "$0")"
for s in 00_setup 10_corrected_corpora 20_rankings 30_controls 40_metastasis 45_comparators 50_downstream 65_corpus_audits 70_fixed_point 60_figures_pdf; do
  echo "=== $s ==="; bash "$d/$s.sh"
done
