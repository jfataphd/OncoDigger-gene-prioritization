#!/bin/bash
# Stage 20: the main 19-cancer ranking, the vocabulary-corrected rankings and their controls, the sensitivity analyses and the proximity runs.
# Longest stage (several hours on 16 logical processors; see docs/REPRODUCTION.md).
. "$(dirname "$0")/env.sh" || exit 1
run a01_run_query_planes scripts/analysis/run_query_planes.py --oncodigger-root "$ONCODIGGER_ROOT" --planes mutation
run a02_vocab_corrected_rescoring scripts/analysis/vocab_corrected_rescoring.py
run a03_vocab_correction_supporting_papers scripts/analysis/vocab_correction_supporting_papers.py
run a04_met_strict scripts/analysis/vocab_corrected_rerun_met_strict.py
run a05_apply_vocab_correction scripts/analysis/apply_vocab_correction_to_ranking.py
# the canonical ranking table (data/canonical/precision_at_k_detail.csv) must be rebuilt from the final ranking BEFORE the scripts below read it
run a06_build_primary_tables scripts/analysis/build_primary_tables.py --ranking data/derived/query_planes/mutation_vocab_corrected.csv
for s in compute_frequency_baseline_with_aliases compute_frequency_baseline_primary_only matched_frequency_baseline intogen_literature_weight_audit true_raw_vs_cleaned_all19 retrieval_vocab_2x2_faithful intogen_full_corpus_presence_check intogen_miss_root_cause; do
  run a10_$s scripts/analysis/$s.py
done
run a20_driver_mutation_plane scripts/analysis/run_query_planes.py --oncodigger-root "$ONCODIGGER_ROOT" --planes driver_mutation
run a21_intogen_parameter_sensitivity scripts/analysis/intogen_parameter_sensitivity.py
run a22_pool_size_jaccard scripts/analysis/pool_size_jaccard.py
# proximity runs use the corrected-collections root
D=data/derived/proximity; PX="$WORK/proximity"; mkdir -p "$PX"
( export ONCODIGGER_ROOT="$WORK/s27corr"
  for w in 0 3 0.25 1 10; do
    mkdir -p "$PX/w$w"
    run a30_prox_w$w scripts/analysis/proximity_boost_run.py $w "$PX/w$w"
    cp "$PX/w$w/mutation.csv" $D/ranking_w$w.csv
    [ "$w" = "3" ] && cp "$PX/w$w/scores_w3.csv" $D/scores_w3.csv
  done
  for w in 0 3; do
    mkdir -p "$PX/rule_w$w"
    run a31_prox_rule_w$w scripts/analysis/proximity_boost_run.py $w "$PX/rule_w$w" rule
    cp "$PX/rule_w$w/mutation.csv" $D/ranking_rule_w$w.csv
  done ) || exit 1
echo "STAGE 20 DONE"
