# Revision after the vocabulary correction (policy v5), 2026-10-17

Revision notes written on 2026-10-17. They build on `docs/REVISION_exact_background_2026-10-16.md` (exact background count). The state before this
revision is not part of this release.

## What was done

* One frozen, versioned vocabulary policy, `scripts/analysis/vocab_policy.py` (version `2026-10-16-v5`), read by every
  cleaned-vocabulary analysis. Decisions come from `data/derived/vocab_rescoring/vocab_policy_v5_table.csv` (122 rows: 69 suppress,
  45 retain, 8 retain-collection-dependent), derived by a rule written before the rebuild and independent of IntOGen.
* Whole-vocabulary screen (`vocab_inventory_screen.py`): 133 alias and 32 symbol candidates; 119 belong to the 602 genes in
  reported rankings and controls. v3: 13 symbols; v4 adds SHE, RHO, CCS, NODAL and the aliases rank, nude, lobe; v5 adds 14 symbols
  and 48 aliases. Final: 7 genes excluded, 31 bare symbols and 51 aliases suppressed at term level (78 of the 82 declared terms
  change the matcher; GC, HR, HCCS and BPIFA4P are not in the base vocabulary).
* Matcher-level preflight test (`test_vocab_policy_v5.py`): 195 checks, 0 failures.
* Full rebuild of every cleaned-vocabulary result on the frozen 2026-09-27 engine with the exact-background hook (invariants
  n <= B <= C checked at every candidate gene; no run stopped): main ranking, MET/REST/FH correction, pooled control, 2022 and 2023
  cutoff rankings, driver-mutation plane, parameter sensitivity, pool-size Jaccard, seven proximity runs, gene-exclusion runs, flat-count
  controls and the retrieval-by-vocabulary factorial. True-raw control arms are unfiltered and were not rerun.
* All dependent tables, figures (16 sources and PNGs), provenance record and manuscript text regenerated from the outputs.

## Before and after (the values that changed)

| Quantity | Before v5 | After v5 |
|---|---|---|
| Any-cancer Precision@10, exact concordance | 178/190, 150/190 | unchanged |
| Top-10 genes | | identical in all 19 cancers (0 of 190 slots changed) |
| Top-50 positions changed (primary ranking) | | 30 of 950 |
| Pooled control (any, exact) | 190/190, 101/190 | unchanged (top 10 identical) |
| Table 2 cell 2, full corpus cleaned (any, exact) | 82, 49 | **108, 62** |
| Table 2 cell 4, BM25 pool cleaned (any) | 174 | 175 |
| Interaction, any-cancer | +44.2 pp (CI 35.8-52.6), positive 19/19 | **+31.1 pp (CI 21.1-41.6), positive 16/19** |
| Interaction, exact | +46.3 pp (CI 38.9-53.7) | **+39.5 pp (CI 31.1-47.9), positive 18/19** |
| Matched counting vs full formula, any-cancer | +2.1 pp, p=0.042 | **+1.6 pp, CI 0.0 to 3.2, p=0.083** |
| Matched counting vs full formula, exact | +2.6 pp, p=0.056 | unchanged |
| Frequency-only control, any / exact | 85 / 58 (44.7% / 30.5%) | **118 / 78 (62.1% / 41.1%)** |
| Mean cumulative precision at k=20, 40, 50 | 87.9, 80.0, 75.9% | 88.7, 81.4, 77.6% |
| Recall at rank 20 and 50 (pooled) | 14.5%, 26.1% | 14.6%, 26.3% |
| Rank AUC mean (range) | 0.700 (0.524-0.846) | 0.706 (0.530-0.886) |
| Spearman significant after Bonferroni | 3 cancers | 2 cancers |
| Unrecovered pairs under the strictest filter | 187 | 191 |
| Hook candidates evaluated; engine B lower | 4,323; 3,544 (82.0%) | 3,955; 3,202 (81.0%) |
| Gene-exclusion test, top-50 slots changed | 203 | 208 |

Unchanged: CancerMine, ChatGPT, DISEASES and PubTator results, the true-raw arms, the cutoff-year headline values, the proximity
null result, the driver-mutation variant (+2.1 and +1.6 points), the pool-size Jaccard (mean adjacent 0.897).

## Interpretation notes

* Table 2 is a retrieval-by-vocabulary factorial with flat mention counting in all four cells. The full relevance formula is
  evaluated only in the separate matched comparison. Cell 4 (175) equals the matched frequency-only baseline (175).
* Cleaning raised the full-corpus condition (82 to 108) but not the BM25-pool condition (174 to 175). This is reported as an
  inference from this experiment, not a general property of BM25.
* The formula's advantage over matched counting is small and not significant; the text says so without claiming equivalence.

## Residual ambiguity (disclosed, no further rebuild)

Targeted check of the 30 gene entries that newly entered a top-50 list: TTR (Colorectal, rank 50; 10 of 10 sampled abstracts were
"time to recurrence"), FH (Breast, 49; 18 of 20 "family history"), APC (Melanoma, 49; 7 of 20 the gene, otherwise
antigen-presenting cells, APC/C, activated protein C, annual percentage change). MKI67, CTLA4, CD34, BAX and BCAR3 were genuine.
Eight collection-dependent gene names (KIT, MET, APC, POLE, DDIT3 via CHOP, MIA, BAD, FAS) and the family aliases (HDAC, VEGFR,
HSP90, PARP, PI3K) are retained and disclosed. The term readings were made by a language-model agent that had seen the
manuscript and are not independent validation. The count of three is not a verified bound on all remaining ambiguous matches.

## Symbol-presence eligibility filter (found 2026-10-18 when the app parity test stopped at Breast rank 47)

`rank_genes` in `run_query_planes.py` skips a gene whose own lower-case symbol has no posting in the cancer's collection. The rule was in the first version of
the script, was not tuned against IntOGen, and was not described in the manuscript. `scripts/analysis/symbol_gate_diagnostic.py` (output
`data/derived/symbol_gate_diagnostic.csv`, reproduces mutation.csv 950/950) shows it removed 18 candidate entries in 10 cancers at engine positions 16 to 54 and none in a
top 10; 16 of them lay inside the 50 places and their replacements entered at ranks 47 to 50 (9 of 16 IntOGen drivers against 0 of the 16 removed). Mean cumulative
precision without the rule would be 0.26, 0.53, 0.92 and 0.95 points lower at k=20, 30, 40 and 50 and unchanged at k=10. It applies to the rankings made by `rank_genes`
(main, driver-mutation, cutoff, pooled control, proximity, relevance ranking of the matched controls) and not to the flat-count rankings of Table 2 or the frequency-only
controls. Now documented in Methods 3.3, Supplementary Methods S1 and Section 7. The analyses were not rerun and the filter was not removed. (An earlier message
in this session said 19 removed entries; the correct count is 18.)

## Fixes found on the way

* `intogen_build_suppfig9_retrieval_vocab_interaction.py` assumed every per-cancer interaction was positive (fixed 0-8 axis);
  it now extends below zero.
* `intogen_build_figure1_workflow.py` takes the vocabulary counts from the policy module.
* `check_provenance.py` expected REST to change Thyroid ranks 41-50; under v5 the lists are identical with and without the redaction.
* `data/provenance/hook_stats_main_ranking.json` replaced with the v5 run's statistics.
* Ordering of `secondary_stats.py` relative to `intogen_precision_at_k_rescoring.py` (secondary_stats must run after).

## Not done, on hold

Open at the time of writing: the release archive (v1.3, Zenodo) and the availability statement, the port of the v5 policy to the live
application (its `gene_quality.py` was at an earlier policy), the term-level PMID export and the short-alias (p53) decision.
