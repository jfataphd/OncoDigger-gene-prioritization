# Build provenance

Every result in the manuscript is computed from **one** corpus build, the **2026-09-27 release** of the
19 cancer-specific PubMed databases (1,640,505 abstracts). This folder records what that build is and ties each
committed output to it.

| File | What it is |
|---|---|
| `corpus_build_2026-09-27.json` | The 19 databases of the build: file name, SHA-256, byte size, document count, plus the pooled index and the engine identity. |
| `raw_reconstruction_2026-09-27.json` | How the raw abstract files needed by the audit scripts were reconstructed and verified against the databases. |
| `outputs_record.csv` | SHA-256 of every committed output (`data/derived`, `data/canonical`, `figures`), the corpus build and the frozen vocabulary policy (`vocabulary_policy`, 2026-10-16-v5) it belongs to, and the script(s) that write it. Outputs with no generator script are labelled `NO GENERATOR IN REPOSITORY`. |

## Build

The 2026-09-27 build is the 2026-05-08 baseline (1,608,566 abstracts) plus an incremental PubMed fetch for
publication dates 2026/05/01 to 2026/09/27, net of records removed during the refresh. The databases are about
27 GB and are not bundled. They are held by the corresponding author. A rebuild from NCBI E-utilities will not be
identical (PubMed records are revised and withdrawn), which is why every database is identified by SHA-256.

## Guard

`scripts/analysis/provenance.py` refuses to run a corpus script on any other build:

* `run_query_planes.py` calls `assert_corpus_build(ONCODIGGER_ROOT)` and stops if any of the 19 databases has a
  different size or document count from the manifest (pass `--allow-nonstandard-build` only for deliberately
  different corpora such as the year-truncated cutoff runs).
* `scripts/analysis/check_provenance.py` checks that the committed outputs are unchanged since recorded, that
  `corpus_stats.csv` equals the manifest, and that independent outputs agree with each other (canonical vs
  corrected ranking, headline counts across seven files, the 1,740-pair partition across three files).
  With `ONCODIGGER_ROOT` set it also verifies the corpus root; add `--full-hash` to recompute every SHA-256.

Regenerate the record after rerunning the analysis, then check it:

```
python -X utf8 scripts/analysis/make_provenance_record.py
python -X utf8 scripts/analysis/check_provenance.py
```

## Engine and environment

* Ranking engine: private live-application repository, commit `4e39fc068102` (2026-05-31). The four engine files
  used are byte-identical to `reference_engine_code/oncodigger_search/`.
* Scripts that read the application's unfiltered gene lookup import `gene_quality.py` policy `2026-09-27-v1`
  (live-application commit `d592370`), vendored in `reference_engine_code/gene_quality_test_engine/` with
  `SHA256SUMS.txt`. The later policy `2026-10-05-v2` described in the manuscript's supplement is a live-application
  policy; the analysis scripts did not use it.
* All engine runs use `PYTHONHASHSEED=0`.

## Known limits

* Three scripts query live services and are date dependent: `cancer_specific_drivers_pubmed_validation.py`
  (NCBI esearch), `pubtator_benchmark_v2_subtype_aggregated.py` and `pubtator_v1_single_parent_reconstruction.py`
  (PubTator). Their committed outputs are frozen snapshots and must not be regenerated for the manuscript.
* `intogen_figure4_headtohead_data.py` needs `data/derived/pubtator_top10.csv`, the original pre-correction
  PubTator list, which is intentionally not in this repository; its output is committed.
* `corpus_stats.csv` has no generator script; it is verified against the manifest by `check_provenance.py`.
