# IntOGen independence audit

Row-level audit behind Supplementary Results S13 of the manuscript: how much of OncoDigger's IntOGen concordance
depends on IntOGen's literature-informed admission rules (Cancer Gene Census and CancerMine), and how the 190
top-10 pairs split into evidence classes A, B and C.

The two scripts were first written as one-off audit copies (2026-10-11) with absolute paths from the author's
machine. They are now repository-relative and runnable:

```
ONCODIGGER_ROOT=<2026-09-27 corpus folder> PYTHONHASHSEED=0 python -X utf8 scripts/analysis/independence_audit/paper1_audit_benchmarks.py
```

(`paper1_audit_benchmarks.py` imports and runs `paper1_audit_main.py`, then recomputes every system's
benchmark under each reference.) Inputs are all in the repository: the IntOGen compendium and unfiltered table
(`data/external/intogen/`), the intogen-plus v2024 whitelist (`data/external/intogen_plus/white_listed.txt`,
`build/datasets/postprocess/white_listed.txt` of `github.com/bbglab/intogen-plus` tag `v2024`), the canonical
ranking and the comparator top-10 files. `ONCODIGGER_ROOT` is needed only for two corpus-volume columns of
`paper1_190_pair_evidence.csv`, which do not affect any class. Outputs are written to
`data/derived/intogen_independence_audit/`.

On the 2026-09-27 build the audit reproduces the independent rescoring script
(`intogen_statistical_support_rescoring.py`) for every system, for example OncoDigger 158/190 any-cancer and
125/190 exact under the statistical-rule reference, CancerMine 152/99, ChatGPT 147/130, pooled control 171/74 and
matched frequency-only 157/128. The two code paths are independent, so agreement is a cross-check.

## Limits

CGC membership is licence-restricted and not released, so CGC status is a lower bound inferred from rule logic.
IntOGen's CGC-calibrated method weights cannot be removed without re-running the IntOGen pipeline on the raw
cohorts, which was not done. The restricted references are sensitivity analyses, not independent references.
The 142 rows removed by the literature filter can be confirmed from `Unfiltered_drivers.tsv`
(`FILTER == "Lack of literature evidence"`); none matches any of the 190 top-10 pairs.
