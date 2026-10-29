# External comparator data provenance

Addresses red-team audit item 18: exact release/version, filename, checksum, and access date
for every external dataset this paper's head-to-head benchmarks depend on. CancerMine, DISEASES,
and PubTator each publish bulk data differently (a versioned, DOI-archived snapshot; a
continuously-updated bulk file with no version number; a continuously-updated bulk file with no
version number, respectively), so "CancerMine" or "PubTator" alone does not by itself identify a
reproducible dataset without the specifics below.

| Resource | Release identity | File | SHA-256 | Access date |
|---|---|---|---|---|
| CancerMine | Zenodo record 7689627 (v50, published 2023-03-01; a fixed, versioned, DOI-archived snapshot — the record itself does not change). Both CancerMine files can be re-downloaded and verified with `scripts/reproduce/fetch_external.py`. | `cancermine_collated.tsv` | `3f5bfc625567439cb0e2e1f96e89d96b149d2d1070457cff1c464563da233638` | 2026-09-26 |
| CancerMine (sentence-level, item 4/7 audits) | same Zenodo record 7689627 | `cancermine_sentences.tsv` | `aa4fd0901d343bdf2987a76a2d60ca8a39d58d73c87333a2b750393403e6a74d` (Zenodo MD5 `defa8f4e07dc8e1400ba363284f562c6`, 83,204,702 bytes; verified 2026-10-18) | 2026-09-30 |
| DISEASES 2.0 | Full bulk "textmining" channel release; Jensen Lab does not publish a version number or fixed release date for this file — it is continuously updated, so the access date below is this dataset's only available identity anchor | `textmining_full_cancer_subset.tsv` (cancer-related subset streamed from the ~1.9GB full release; not committed, gitignored) | `ddd3cb6e059ee0df31d2d4f9776864f1220869ed6375f4233161db361e2ebef9` | 2026-09-26 |
| PubTator 3.0 | `relation2pubtator3.gz` bulk relation-extraction release; NCBI does not publish a version number or fixed release date for this file — it is continuously updated, so the access date below is this dataset's only available identity anchor | `relation_cancer_subset.tsv` (filtered to the 19 cancers' MeSH IDs; not committed, gitignored) | `ccd1162be2b1d61344b3737abdb89ce773fa66dc0c6bb55e2f20f6e0efed66c9` | 2026-09-26 |
| PubTator 3.0 gene-ID mapping | NCBI `Homo_sapiens.gene_info.gz` (continuously updated, no version number) | `gene_id_to_symbol.csv` (small derived Entrez-to-symbol lookup, committed) | `dcce32f8432f8f33fb20bedff9bc5c99f75e419e08927494b091004bf3cc82e7` | 2026-09-26 |
| English common-word list (20,000 words) | github.com/first20hours/google-10000-english, `20k.txt`; not shipped with this release, `scripts/analysis/vocab_inventory_screen.py` fetches it when absent and refuses a copy with a different checksum | `data/reference/english_common_words_20k.txt` | `4ed6e5336d7760d281f7e72df31827da880c861363e820d8c65666b0f10d9ac0` (155,418 bytes) | 2026-10-29 (byte-identical to the copy the paper's vocabulary screen was run with) |

**Temporal-fairness note (item 18's own concern):** because DISEASES and PubTator publish
continuously-updated bulk files with no version identifier, "CancerMine/DISEASES/PubTator" as of
2026-09-26 are not necessarily the same underlying snapshot a future re-download would produce.
This is a genuine, unavoidable limitation of benchmarking against these two resources specifically
(unlike CancerMine's fixed Zenodo record), stated here rather than left implicit; the checksums
above are the definitive record of exactly what data these results are based on regardless of what
either resource looks like at some future access date.

**Archiving plan for the raw files themselves.** The three gitignored files above (cancer-relevant
subsets, not the full upstream releases) are 83MB (CancerMine sentences), 269MB (DISEASES) and
49MB (PubTator) -- the DISEASES file alone exceeds GitHub's 100MB per-file limit without Git LFS,
so none of the three can simply be committed to this repository. The checksums above let anyone
verify a re-download matches exactly what these results are based on; for the raw bytes themselves,
the plan is to upload all three alongside the frozen manuscript release when this repository is
archived to Zenodo (the same standing pre-submission step already gating public release), rather
than introduce Git LFS or separate hosting infrastructure just for this.
