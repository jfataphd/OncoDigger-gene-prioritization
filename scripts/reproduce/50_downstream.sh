#!/bin/bash
# Stage 50: everything computed from the committed ranking and comparator files (no corpus needed): tables, statistics, bootstraps, figure data,
# comparator scoring, sensitivity summaries, the headline verification and the manuscript number check. About 15 minutes.
# Expected result: every committed file regenerates byte-identical (check with `git status`: no modified files).
export REPRO_NO_CORPORA=1
. "$(dirname "$0")/env.sh" || exit 1
unset ONCODIGGER_ROOT BGFIX_ENABLE
r() { run "$@"; }
"$PYTHON" -X utf8 scripts/reproduce/fetch_external.py || exit 1   # CancerMine files (fixed Zenodo record), checksum-verified
r d01_build_primary_tables scripts/analysis/build_primary_tables.py --ranking data/derived/query_planes/mutation_vocab_corrected.csv
r d02_intogen_full_rescoring scripts/analysis/intogen_full_rescoring.py
r d03_precision_at_k_rescoring scripts/analysis/intogen_precision_at_k_rescoring.py
r d04_secondary_stats scripts/analysis/secondary_stats.py
r d05_statistical_support scripts/analysis/intogen_statistical_support_rescoring.py
r d06_interaction_all19 scripts/analysis/interaction_all19.py
r d07_interaction_bootstrap scripts/analysis/interaction_all19_bootstrap.py
r d08_figure3_data scripts/analysis/intogen_figure3_data.py
r d09_corpus_specific_genes scripts/analysis/corpus_specific_genes_full_list.py
r d10_cutoff_rescoring scripts/analysis/cutoff_matched_rescoring.py
r d11_gene_exclusion_compare scripts/analysis/gene_exclusion_readmission_test.py compare
r d12_proximity_summary scripts/analysis/proximity_summary.py
r d13_proximity_vs_matched scripts/analysis/proximity_vs_matched_counting.py
r d14_mapping_code_sensitivity scripts/analysis/mapping_code_sensitivity.py
r d15_driver_mutation_sensitivity scripts/analysis/driver_mutation_sensitivity.py
r d16_figure4_headtohead scripts/analysis/intogen_figure4_headtohead_data.py
r d17_precision_bootstrap_ci scripts/analysis/intogen_precision_bootstrap_ci.py
r d18_cleaned_recall scripts/analysis/intogen_cleaned_recall.py
r d19_gene_breadth scripts/analysis/intogen_gene_breadth.py
r d20_confidence_breakdown scripts/analysis/intogen_confidence_breakdown.py
r d21_cancermine_dedup scripts/analysis/intogen_cancermine_dedup_exact_concordance_check.py
r d23_pooled_exact scripts/analysis/pooled_control_exact_concordance.py
r d24_ablation_panel scripts/analysis/intogen_ablation_panel.py
r d25_metastasis_summary scripts/analysis/metastasis_exclusion_summary.py
r d26_metastasis_headline scripts/analysis/metastasis_exclusion_headline.py
r d27_rank_aware_metrics scripts/analysis/rank_aware_metrics.py
r d28_verify_headline scripts/analysis/verify_headline_results.py
r d28b_export_repository_tables scripts/analysis/export_repository_tables_rt10_rt12.py
r d29_check_text_numbers scripts/manuscript/check_text_numbers.py
echo "STAGE 50 DONE"
