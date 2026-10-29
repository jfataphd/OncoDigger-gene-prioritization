# Vendored engine source (reproducibility pin)

Addresses red-team audit item 15: the manuscript's Methods described the "Primary Evidence"
weighting function and the BM25 retrieval implementation as living in a separate application
repository (`github.com/jfataphd/OncoDigger`), which continues to be actively developed. Citing
that repo without an exact commit pin means a future reader cloning it would get a different,
evolved version of this code, not the one that actually produced this paper's numbers.

These ten files (all of `src/oncodigger/search/`) are an exact, byte-for-byte copy (as committed, with LF line endings; a Windows checkout with core.autocrlf converts them to CRLF, which does not affect the code) of the engine source this paper's entire
analysis pipeline (`scripts/analysis/run_query_planes.py` and everything downstream of it)
actually imports and calls, at the exact commit the manuscript's own frozen release depends on:

- **Source repository:** `github.com/jfataphd/OncoDigger`
- **Exact commit:** `4e39fc0681024dc4cef5b28978386f40878bcb4e`
- **Commit date:** 2026-05-31
- **Commit message:** "Add alias substitution for ambiguous gene symbols and gene validation UI"
- **Verified:** working tree at this commit has zero uncommitted local modifications (`git status`
  clean) as of the date these files were copied into this repository, 2026-09-30.

Files (original path `src/oncodigger/search/` in the source repository):

| File | Role in this paper's pipeline |
|---|---|
| `lexical.py` | BM25 retrieval (k1=1.5, b=0.75) and metadata multipliers (title/MeSH/recency/review) |
| `discovery_modes.py` | The "Primary Evidence" content-type weighting function (Methods 2.1) |
| `discovery.py` | Entity aggregation and gene enrichment scoring |
| `lexicons.py` | HGNC gene vocabulary loading, including the `useful_gene_term` common-word filter |
| `association.py` | `enrichment_score` and the engine's original, approximate `entity_background_papers` (superseded for every reported result by the exact count, see below) |
| `article_types.py`, `embeddings.py`, `semantic.py`, `schema.py`, `__init__.py` | Imported by `discovery.py` (not used by the mutation-focused ranking, but required for the package to import) |

These are reference copies for reproducibility only; they are not re-imported by anything in
this repository (the analysis scripts still import live from `--oncodigger-root`, defaulting to
the frozen local working copy at the commit above). If that working copy is ever lost or
modified, these files are the fallback source of truth for exactly what code produced every
number in this paper.

Also included here, extracted from the same exact commit for the same reason: the vocabulary
files `lexicons.py` loads. `Human Genes.csv` is the HGNC gene vocabulary (HGNC data is CC0).
`load_lexicons()` loads four files unconditionally, and only the gene vocabulary is used by the
mutation-focused ranking this paper reports.

The other three lists, `kegg_drug_list_lowercase.csv`, `kegg_compounds_lowercase.csv` and
`phytochemicals.csv`, are not redistributed: KEGG's terms do not grant redistribution, and the origin
of the phytochemical list is not recorded. They are shipped as header-only placeholder files (headers
`Drug-ID,drugs`, `Compound ID,compound` and `phyto`; the real files hold one lowercase name per row,
synonyms separated by `; `), so that `load_lexicons()` still starts. This does not change what the
paper reports: the exact-background ranking of Endometrial Cancer, Myeloma and Bladder Cancer (150
rank positions), its exact-background statistics, and both frequency baselines came out identical
with the placeholders and with the real lists (tested 2026-10-29). The lists matter only for
rebuilding a collection from its PMIDs with the original tokenisation of multi-word drug, compound
and phytochemical names (`scripts/reproduce/rebuild_corpus_from_pmids.py`); an exact rebuild needs the
original lists. The KEGG lists can be obtained from KEGG under its own terms and placed here; the
phytochemical list is not available from this release.

Re-running `run_query_planes.py` end to end additionally requires the 19 cancer-specific
SQLite corpora it builds its ranking from (Methods 2.1) — not included here, since they total
over 1.6 million PubMed abstracts; they are reproducible from the freely available NCBI
E-utilities API (exact queries: Supplementary Table S2). Every result built on top of the
frozen ranking this script already produced (`data/canonical/precision_at_k_detail.csv`) —
the IntOGen scoring, all comparator benchmarks, all ablations and controls — is fully
reproducible directly from the files in this repository, independent of the raw corpora.

## Exact background count (added 2026-10-16)

Every result in the manuscript uses the exact background count of `docs/BACKGROUND_COUNT_SPEC.md`, not `association.entity_background_papers`
in the file above. The count is supplied at run time by `exploratory/background_rebuild/hook/` (`bgfix_hook.py`), which replaces the
background-count step of `discovery.associated_entity_tables` in memory when `BGFIX_ENABLE=1`; the files here are not edited. At every
candidate gene of every run the hook checks that the pool support equals its own replica count, that support <= background <= collection size,
and that 0 < enrichment <= collection size / pool size, and it stops the run on any violation. Pre-v5 per-collection count files from the pilot are superseded; see `exploratory/background_rebuild/README_superseded.md`. The package in this folder imports on its own
(`sys.path.insert(0, "reference_engine_code"); import oncodigger_search.discovery`); the commit above is the engine commit, so the vendored
files and the corpus-root copy used for the runs are identical (checked by SHA-256 with carriage returns removed).
