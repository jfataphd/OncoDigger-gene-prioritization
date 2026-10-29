#!/bin/bash
# Stage 30: pooled-corpus control, publication-year cutoff runs (2022, 2023) and the gene-exclusion robustness tests.
. "$(dirname "$0")/env.sh" || exit 1
run b00_pooled_index scripts/analysis/build_pooled_index.py
run b01_pooled_corpus_control scripts/analysis/pooled_corpus_control.py
for y in 2022 2023; do
  if [ ! -f "$WORK/trunc/y$y/build_report.csv" ]; then
    run b01_truncate_$y scripts/analysis/build_date_truncated_corpora.py --src-root "$ONCODIGGER_ROOT" --out-root "$WORK/trunc/y$y" --year-max $y
  fi
  "$PYTHON" -X utf8 scripts/reproduce/assemble_roots.py --engine-into "$WORK/trunc/y$y" || exit 1
  mkdir -p "$WORK/cut$y"
  run b02_cut_rank_$y scripts/analysis/run_query_planes.py --oncodigger-root "$WORK/trunc/y$y" --allow-nonstandard-build --planes mutation --out-dir "$WORK/cut$y"
  cp "$WORK/cut$y/mutation.csv" data/derived/cutoff_matched_rankings/mutation_corpus_to_$y.csv
done
run b03_cutoff_vocab_correction scripts/analysis/cutoff_matched_vocab_correction.py --trunc-root-2022 "$WORK/trunc/y2022" --trunc-root-2023 "$WORK/trunc/y2023" --work-dir "$WORK/cutcorr"
for v in A B C; do run b10_gene_exclusion_$v scripts/analysis/gene_exclusion_readmission_test.py run $v; done
echo "STAGE 30 DONE"
