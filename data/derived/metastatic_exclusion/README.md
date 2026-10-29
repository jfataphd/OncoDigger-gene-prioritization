# Post hoc sensitivity analysis: records indexed under the cancer term's /secondary MeSH subheading (Supplementary Results of the manuscript, Section 7)

Question (from review): the collections are built with `(<MeSH term>[mh]) AND english[la] AND hasabstract`, and `[mh]` includes records indexed with the /secondary subheading,
among them studies of metastases to the named site. How sensitive are the rankings and the headline concordance to removing those records?

What was done, for all 19 collections, with the frozen 2026-09-27 corpora, policy v5, the symbol-presence filter, unchanged BM25 settings and formula:
1. `scripts/analysis/fetch_secondary_pmids.py` lists the PMIDs indexed as `"<term>/secondary"[mh]` (PubMed E-utilities run on 2026-10-17/18, yearly windows) and intersects them with the frozen
   collection: `secondary_pmids_<collection>.csv`, `secondary_summary.csv` (percentages are of the frozen collection sizes).
2. `scripts/analysis/build_metastasis_excluded_corpora.py` builds copies of the collections without those records and recomputes N, mean length, IDF and postings
   (formulas of `build_date_truncated_corpora.py`): `build_report_metastasis_excluded*.csv`. The recomputed tables were verified against the data
   (`build_date_truncated_corpora.py --verify`): `integrity_check_*.log`, all passed.
3. `scripts/analysis/run_query_planes.py --oncodigger-root <reduced root> --allow-nonstandard-build --planes mutation --cancers ...` with the exact-background hook, in two batches
   (Liver, Brain, Lung, Kidney, Melanoma into `ranking/`; the other 14 into `ranking_other14/`): `run_times_batch*.txt`.
4. `scripts/analysis/metastasis_exclusion_summary.py` and `metastasis_exclusion_headline.py` compare with the main ranking (`data/derived/query_planes/mutation.csv`):
   `comparison_per_cancer.csv`, `gene_movements_top10.csv`, `liver_selected_genes.csv`, `comparison_overall.csv`, `headline_main_vs_excluded.csv`.

Result: 6 of 190 top-10 positions changed in five cancers; any-cancer Precision@10 178/190 -> 177/190, exact concordance 150/190 -> 150/190. Leukemia, Lymphoma and Myeloma (0, 2 and 70 records removed)
reproduce their main top-50 lists exactly. The pooled control and the external comparators were not recomputed.

Limits: /secondary is an imperfect proxy (it removes some unrelated and some relevant primary-cancer papers and misses metastasis papers indexed otherwise); for Melanoma the records concern
metastatic melanoma itself. This is a sensitivity analysis, not a primary-tumour classification or an estimate of true contamination. The main analysis, the corpora and all primary results are unchanged.
