# Gene-recognition audit of the top-10 evidence (Supplementary Results S16, Supplementary Table S13)

Post hoc audit of which gene each matched term credits in the supporting abstracts of the 190 top-10 gene-cancer pairs, and read-only
counterfactual rescorings under alias removals. It was designed after the rankings were frozen, used one language-model reader that was
not blinded, and has no independent human validation (see Supplementary Results S16). It does not assess missed mentions, genes ranked
11 to 50, or whether an abstract reports a mutation in the credited gene.

No PubMed abstract text is redistributed here. The per-credit files give PubMed identifiers and titles only.

## Files

| Folder | Contents |
|---|---|
| `per_credit/<cancer>_qc_final.csv` | One row per gene-abstract credit that was read (8,489 rows over 19 cancers). Columns: `cancer, gene, pmid, sampling_method` (random, all-matched or targeted), `triggering_alias`, `identity` (C correct, X incorrect, U uncertain), `evidence_draft`, `target_context`, `rationale`, `pool_rank`, `title`, `evidence_final`. |
| `per_credit/<cancer>_verdict_overrides.csv` | Reader verdicts that replaced the rule-based draft for some credits (where present). |
| `summaries/ALL19_framework_pair_rates.csv` | Per gene-cancer pair: credits, random-sample size and incorrect count with Wilson interval, confirmed incorrect (lower bound) and uncertain counts. Sums: 22,071 credits, 7,808 random or all-matched, 89 incorrect, 213 confirmed by rule. |
| `summaries/ALL19_poolwide_rule_counts.csv`, `ALL19_affected_gene_rates.csv`, `ALL19_error_mechanisms.csv` | Whole-pool counts for each recurring error mechanism (alias rules checked against the abstracts that were read). |
| `summaries/ALL19_evidence_random_or_all.csv`, `ALL19_evidence_targeted_only.csv`, `ALL19_cancer_summary.csv`, `ALL19_suspected_false_matches.csv` | Per-pair evidence-class counts, per-cancer summary and the credits judged incorrect or uncertain. |
| `summaries/<cancer>_qc_gene_summary.csv` | Per-gene read, random, incorrect and uncertain counts for each cancer. |
| `counterfactuals/cf_<cancer>_<variant>.json` | Top-50 ranking (rank, gene, relevance, supporting papers, exact background count) after removing one alias mapping from the frozen vocabulary, with the unmodified baseline `cf_<cancer>_baseline.json` from the same harness. |
| `counterfactuals/CF_consolidated*.csv` | Before and after rank, supporting count, exact background, relevance and top-10 replacements for the 17 interventions (Supplementary Table S13). |
| `counterfactuals/CF_baseline_verification.csv` | Each baseline compared with the frozen v5 top 50: gene order for all 50 positions, and supporting count, exact background and relevance for the IntOGen driver genes among them (the genes listed in `data/derived/intogen_rescoring/intogen_literature_weight_audit.csv`, 13 to 31 per cancer), plus IntOGen 178/190 and 150/190. All equal except Lung MET, which that file records before the MET correction (67 supporting papers; 58 after). |
| `counterfactuals/CF_all_position_changes.csv` | Every changed rank position in every counterfactual run. |

## What can and cannot be re-run

The scripts are in `scripts/analysis/gene_recognition_audit/`. They are the as-run copies; only absolute paths were replaced (the repository
root comes from `ONCODIGGER_PAPER_REPO` or from the script location). They were run in a work folder that held the scripts, the per-cancer
reading sets and an assembled corpus root (`root/`); reading sets contain abstract text and are not distributed.

* Sampling (`build_set.py`, seed 42), the rule-based draft classes (`classify.py`), the merge of reader verdicts (`finalize.py`) and the
  pooled summaries (`synthesis.py`, `apply_framework.py`, `poolwide_rules.py`) are deterministic given the reading sets.
* The identity verdicts themselves are the reader's labels and cannot be recomputed. The per-credit files are the record of them.
* The counterfactual harness (`counterfactual.py`, `counterfactual_breast.py`, `score_cf.py`) re-runs retrieval and gene
  recognition on the frozen 2026-09-27 collections with the exact background count, so it needs those databases.
* `analyze_cf.py` does not need the databases. Run from a folder holding the committed `cf_*.json` files (with `ONCODIGGER_PAPER_REPO` set to the
  repository if the script is run from elsewhere), it rewrites `CF_consolidated.csv`, `CF_all_position_changes.csv` and `CF_baseline_verification.csv`;
  on the committed JSON files the first two come out byte-identical.
