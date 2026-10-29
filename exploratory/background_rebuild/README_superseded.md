# exploratory/background_rebuild: what is final and what is superseded

**Final (used by every v5 result):** `hook/` (`bgfix_hook.py`, `bgfix_worker.py`, `sitecustomize.py`). It computes the exact background count at run time
when `BGFIX_ENABLE=1`; `scripts/reproduce/env.sh` puts it on `PYTHONPATH`. Its per-run check statistics are in `data/provenance/hook_stats_main_ranking.json`.
Method: `docs/BACKGROUND_COUNT_SPEC.md`.

**Historical pilot files kept in this release (pre-v5, not inputs to any result):**
`PILOT_REPORT.md`, `bg_matcher.py`, `exact_background.py`, `test_bg_matcher.py`. They document the Myeloma and Lung pilot that tested the matcher
(224 of 224 and 135 of 135 supporting counts reproduced) before the exact count was adopted. They ran under the vocabulary and matching policy that
preceded policy v5 (matching-policy hash `d2ab7cd4`; frozen policy v5 is `2026-10-16-v5`, hash `30ec0011`).

**Superseded and not distributed in this release:** the per-collection exact counts `b_exact_<collection>.csv` and `.stamp.json` (19 collections),
the run logs, `ranking_bgfix.csv`, `bgfix_percancer.csv`, `bgfix_report.md`, and the scripts that wrote or read them (`rerun_main_ranking.py`,
`bgfix_report.py`, `review_myeloma_genes.py`). They were built on 2026-10-16, before vocabulary policy v5 (2026-10-17 commit). Their top-50 lists
differ from the final v5 lists in 15 of the 19 cancers (genes whose aliases v5 suppresses, for example AMACR via the alias "race" in Kidney Cancer:
23 supporting and 630 background papers in the superseded file, 18 and 213 under v5). No reproduction stage reads them. They remain in the Git
history of the development repository (tag `pre-v5-background-artifacts`).

The final v5 rankings were checked against a fresh unmodified rerun of all 19 collections: top-50 gene order identical, and rank, supporting count,
background count and relevance identical for every gene listed in `data/derived/intogen_rescoring/intogen_literature_weight_audit.csv` except Lung
MET, which that file records before the MET correction (67 supporting papers; 58 after).
