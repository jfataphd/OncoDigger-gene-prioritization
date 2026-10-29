# Relevance-formula study (exploratory, 2026-10-15)

Question: what do the terms of R = W * ln(1+E) * ln(1+S) contribute, and how sensitive are the top-10 lists to the form of the score?
Not part of the manuscript. IntOGen agreement is a sensitivity measure here, not a selection criterion.

Files
- `extract_features.py` (needs the engine and the corrected 2026-09-27 databases; run once) -> `features_raw.csv` (3,892 candidate genes with 2 or more supporting papers, 19 cancers).
  Check: recomputed order reproduces `data/derived/query_planes/mutation_vocab_corrected.csv` exactly (0 of 950 top-50 slots differ).
- `formula_variants.py` -> `variant_summary.csv`, `variant_top10.csv`. 19 scoring variants on the same pools and candidate genes.
  Check: the flat-count variant gives 174/190 and 145/190, the paper's matched counting control.
- `factor_diagnostics.py` -> `factor_diagnostics.csv`. Which factor orders the genes.

Findings (19 cancers, one query, one build)
1. W (rank-weighted sum) carries about 71 to 79% of the variance of ln R among the top candidates, ln(1+S) 16 to 24%, ln(1+E) about 5%.
   Enrichment hardly orders the leading genes (Spearman with R about 0.12).
2. W and S are nearly redundant (Spearman 0.89 to 0.94), so support is counted twice in practice.
3. No variant is significantly better than the baseline on exact concordance (paired cancer bootstrap). W alone and S alone give 145/190 against 142/190.
4. Variants that fail: enrichment alone (48/190 exact), dropping W (129), 1/r^2 decay (111), 1/r decay (132), and the de-duplicated (W/S) form (110).
   So rank weighting must not be too steep and support cannot be removed, while enrichment is only a modifier.
5. Alternative enrichment forms (Fisher, G statistic, ln lift, IDF style) land at 143 to 145/190 and are slightly more cancer-differentiated (pairwise Jaccard 0.19 to 0.20 against 0.187).
Open: the engine's background count B is wrong for genes matched only as phrases (dashed symbols) and uses the smallest alias count; this is not yet corrected in these variants.

Background-count finding (`background_union.py`, `background_effect.py`; exploratory, approximate)
- The engine's background count B is the smallest postings count among the matched single-word terms, 0 for phrase-only genes (fallback to the gene's own support).
  Against B_union (papers containing any active single-word alias in the postings, phrase-only genes by text scan), the engine's B is smaller for 1,553 of 3,892 candidate genes (B_union more than 2x B for 688).
- Rescoring with B_union (same W and S): any-cancer 177 -> 178 and exact 142 -> 144 of 190; 10 top-10 slots change in 8 cancers; 52 of 950 top-50 slots change.
- CORRECTION (2026-10-15): the first run of this script reported 179 and 148. That was an error in the script, which looked up only the pre-remap entity id and gave TP53 and ERBB2 a background of 0
  (then clamped to their own support, maximum enrichment). Fixed by looking up both ids. The 148 figure should not be used.
- B_union is still approximate (it uses postings tokens, which keep hyphens, not the engine's own matcher). The planned rebuild computes B with the same matcher as the pool.
- Looked at together with IntOGen agreement, so a sensitivity result. Not in the manuscript.
