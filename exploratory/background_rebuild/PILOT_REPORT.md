> **Historical, pre-v5.** This pilot ran before vocabulary policy v5. `b_exact_Myeloma.csv` and `b_exact_Lung_Cancer.csv` named below are superseded and are not distributed (see `README_superseded.md`). No v5 result depends on this file.

# Pilot of the exact background count (Myeloma, Lung Cancer), 2026-10-15

Script: `exact_background.py` (engine matcher applied to every paper). Outputs: `b_exact_Myeloma.csv`, `b_exact_Lung_Cancer.csv`. Spec: `docs/BACKGROUND_COUNT_SPEC.md`.
Corpus: corrected 2026-09-27 databases. Run time: Myeloma 26,591 papers 198 s; Lung 178,908 papers 1,435 s (about 125 papers per second, one process).
Estimate for all 19 collections (1.64 million papers): about 3.5 hours in one process.

| Test | Myeloma | Lung |
|---|---|---|
| 1. Replica reproduces engine S (candidates) | 224 of 224 | 135 of 135 |
| 2. S <= B and B <= collection size | 0 violations | 0 violations |
| 6. 0 < E <= C/N for every gene with S > 0 | holds | holds |
| Engine B equals exact B | 54 of 224 | 10 of 135 |
| Engine B smaller than exact B | 170 | 125 |
| Engine B larger than exact B | 0 | 0 |
| Engine B more than 2x smaller | 68 | 43 |
| Engine B below S (floor binds) | 19 | 1 |

Examples (S, engine B, exact B): Myeloma SDC1 34, 5, 861. Lung NKX2-1 14, 14, 964. Lung ERBB2 68, 639, 2,046.

Effect on the ranking (shown to audit the gene changes; not used to judge the correction):
- Myeloma top 10: CCND1 and SDC1 leave, IRF4 and TRAF3 enter; 5 of 50 top-50 genes change.
- Lung top 10: same ten genes, order changes (CD274 moves from 3rd to 8th).

Not yet done: tests 3 to 5 and 7 on a hand-built corpus, version-stamped cache, supporting-abstract audit of changed genes, timing for the other 17 collections.
