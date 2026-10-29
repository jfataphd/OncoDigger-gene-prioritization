#!/bin/bash
# Stage 65: corpus-based audits and sensitivity checks reported in the manuscript's supplement (evidence-audit worksheet, symbol-gate diagnostic,
# vocabulary audit tables, true-raw flat counting, formula x vocabulary interaction, stop-word alias check, rank vs enrichment ablation).
# Needs the corpus databases. Steps that read case-preserved raw abstracts run only when ONCODIGGER_RAW is set; otherwise they are reported as SKIPPED.
. "$(dirname "$0")/env.sh" || exit 1
for s in evidence_audit_worksheet symbol_gate_diagnostic vocab_correction_v4_evidence vocab_inventory_screen true_raw_flat_count_2x2 formula_interaction_dirty_vocab stopword_alias_sensitivity_check rank_vs_enrichment_ablation; do
  run v_$s scripts/analysis/$s.py
done
run v_independence_audit scripts/analysis/independence_audit/paper1_audit_benchmarks.py   # also runs paper1_audit_main.py; needs the corpus root for two volume columns
if [ -n "$ONCODIGGER_RAW" ] && [ -d "$ONCODIGGER_ROOT/data/pubmed_raw" ]; then
  run v_title_case_blind_spot_check scripts/analysis/title_case_blind_spot_check.py
  run v_vocab_collision_risk_screen scripts/analysis/vocab_collision_risk_screen.py
else
  echo "SKIPPED title_case_blind_spot_check, vocab_collision_risk_screen: need the raw abstracts (set ONCODIGGER_RAW before stage 00)" | tee -a "$LOGS/progress.txt"
fi
echo "STAGE 65 DONE"
