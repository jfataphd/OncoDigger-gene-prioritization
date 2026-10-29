# Reproducing the OncoDigger paper 1 results

Every number, table and figure in `manuscript/OncoDigger_manuscript1.1.md` is produced by the scripts in `scripts/` through the stages in `scripts/reproduce/`.
Run on the final commit these stages reproduce the committed outputs (see "What was tested" for the exact record).

## What you need

| Input | Where it comes from |
|---|---|
| This repository | the release archive |
| The 19 corpus databases of the 2026-09-27 build (about 27 GB) | **not in the repository**; SHA-256 of each is in `data/provenance/corpus_build_2026-09-27.json`; PMID manifests are in the repository. Point `ONCODIGGER_CORPORA` at the folder holding the 19 `<collection>_cancer_inclusive_all/lexical_index.db` folders |
| DISEASES 2.0 and PubTator 3.0 cancer subsets (`data/external/diseases/`, `data/external/pubtator/`) | the Zenodo archive of this release (257 MB and 47 MB; not in the GitHub repository, whose file limit is 100 MB); SHA-256 in `data/external/PROVENANCE.md`, and stage 45 refuses to run on any other file |
| CancerMine collated and sentence tables | fetched and checksum-verified by stage 45 from the fixed Zenodo record |
| English common-word list (20,000 words), for the vocabulary screen in stage 65 | not in the release; `vocab_inventory_screen.py` fetches it from github.com/first20hours/google-10000-english and checks its SHA-256 (`data/external/PROVENANCE.md`) |
| Case-preserved raw PubMed abstracts (about 6.5 GB), optional | not in the repository. Needed only for two audit scripts (`title_case_blind_spot_check`, `vocab_collision_risk_screen`) and to recompute the MET/REST/FH decision files; those decisions are committed (`data/derived/intogen_rescoring/vocab_correction_pmids/`), so every other step runs without them. Set `ONCODIGGER_RAW` before stage 00 to include them |

Machine used: Windows 11, 16 logical processors, 63 GB RAM, Python 3.12, Git Bash. `BGFIX_WORKERS` (default 8) sets the processes used for the exact background count.

## Commands

```
export ONCODIGGER_CORPORA=/path/to/pubmed_processed
bash scripts/reproduce/run_all.sh          # all stages in order
```

`env.sh` documents the variables (`REPRO_WORK`, `BGFIX_WORKERS`, `REPRO_RESUME=1` to skip steps already finished). Every step writes a log and a line `name exit=<code> seconds=<n>` to `$REPRO_WORK/logs/progress.txt`; a stage stops at its first failure.

| Stage | What it does | Needs corpora | Measured time |
|---|---|---|---|
| 00 setup | assembles the frozen corpus root and checks corpus hashes | yes | seconds |
| 10 corrected corpora | applies the frozen MET/REST/FH correction to Lung, Kidney, Thyroid | yes | seconds |
| 20 rankings | main and driver-plane rankings, cutoff rankings, frequency baselines | yes | about 15 h wall clock |
| 30 controls | pooled control, 2x2 retrieval x vocabulary, gene-exclusion tests | yes | about 4.8 h |
| 40 metastasis | metastatic-literature exclusion test: corpora without /secondary-indexed records, and their rankings | yes | about 1.3 h |
| 45 comparators | CancerMine, DISEASES, PubTator, ChatGPT top-10 lists and IntOGen scoring | no | under 1 min |
| 50 downstream | every table, statistic and number check built from the rankings | no | about 1 min |
| 65 corpus audits | evidence worksheet, symbol-gate diagnostic, vocabulary audits, true-raw count, formula x vocabulary interaction, stop-word alias check, ablation, independence audit | yes | about 1.9 h |
| 70 fixed point | reruns stage 50 and requires that nothing changes | no | about 1 min |
| 60 figures and PDF | figures, PDF, citation check, reference-existence check, number check, provenance check | no | under 1 min |

Stages 00 to 40 together took about 22.6 hours of step time in the recorded run (sum of the per-step `seconds`; stage 20 also includes one failed first attempt that is not counted). The times are from one machine and depend on `BGFIX_WORKERS`.

## Rebuilding a corpus without the databases

The databases hold PubMed abstracts and are not redistributed. Each collection's PMID list is in `data/provenance/pmid_manifests/`, and one script rebuilds a collection from that list alone:

```
python -X utf8 scripts/reproduce/rebuild_corpus_from_pmids.py myeloma --out <folder> [--api-key KEY] [--compare <stored lexical_index.db>]
```

`bash scripts/reproduce/rebuild_all_corpora.sh <folder> [API key]` does this for all 19 collections (not run by the authors; resumable). The script downloads exactly those records from NCBI (efetch by PMID, so the paper set is fixed although PubMed has grown), cleans and indexes them with the application's steps, and writes a `lexical_index.db` in the schema the ranking scripts read. `--compare` (for anyone who holds the stored database) reports how many records differ.

Tested on Myeloma (26,591 PMIDs, 247 s without an API key): 26,590 records came back (one is no longer retrievable with an abstract); 98.70% of them are identical to the stored database in title, abstract, MeSH list, year and length; there were no abstract or year differences, 316 MeSH lists, 3 titles and 28 document lengths differed. Running `run_query_planes.py` on the rebuilt database gave the same 50-gene Myeloma ranking as the committed one (50 of 50 positions). That test used the original KEGG drug and compound lists and the phytochemical list, none of which this release redistributes (header-only placeholders are shipped; see `reference_engine_code/README.md`); a rebuild with the placeholders tokenises multi-word drug, compound and phytochemical names differently from the original databases, so it will not reproduce them exactly. The other 18 collections have not been rebuilt this way. PubMed records are revised after the fact, so a rebuild can differ slightly from the databases used in the paper; the other stages that need the corpora (the exact background counts and everything built on them) have not been rerun on rebuilt databases.

## Final v5 outputs and superseded intermediates

Every committed file under `data/derived`, `data/canonical` and `figures/png` is tied to the 2026-09-27 corpus build and to frozen vocabulary policy
`2026-10-16-v5`; both are recorded per file in `data/provenance/outputs_record.csv` (columns `corpus_build`, `vocabulary_policy`). The exact background
count is computed live by `exploratory/background_rebuild/hook/` (stage scripts set `BGFIX_ENABLE=1`). The earlier per-collection files
(`b_exact_*`, `ranking_bgfix.csv` and related) were built before policy v5, are superseded, are not distributed with this release and are not read by any
stage; `exploratory/background_rebuild/README_superseded.md` lists what is final, what is kept as a historical pilot, and what was removed.

## Not rerun by the stages (live services or private data)

`pubtator_benchmark_v2_subtype_aggregated`, `pubtator_v1_single_parent_reconstruction`, `cancer_specific_drivers_pubmed_validation` and `fetch_secondary_pmids` query live services whose answers change; their outputs are committed and checksummed (`data/provenance/outputs_record.csv`). The PowerPoint-to-PNG export of the figures needs desktop PowerPoint; the PNGs are committed and recorded in the provenance file. Scripts under `exploratory/` are not part of any reported result unless a stage names them.

## Random seeds and bootstrap procedure

Every confidence interval in the paper is a case bootstrap over the 19 cancers (the cancer type is the unit), with 2,000 resamples and seed 42.

* **Single-system and paired intervals** (`scripts/analysis/intogen_full_rescoring.py`, `intogen_precision_bootstrap_ci.py`, `intogen_statistical_support_rescoring.py`): `numpy.random.default_rng(42)` is created afresh for each statistic and draws one vectorised index array `integers(0, 19, size=(2000, 19))`; the per-cancer values are taken in the study-cancer order of Supplementary Table S5 (`STUDY_ORDER`), the statistic is the mean over each resample (for a paired difference, the mean of the per-cancer differences), and the interval is the 2.5th and 97.5th percentile (numpy default linear interpolation). A per-replicate loop with the same seed draws a different sequence and gives slightly different bounds.
* **Retrieval x vocabulary interaction** (`scripts/analysis/interaction_all19_bootstrap.py`): Python `random.Random(42)`, one `rng.choice(cancers)` per cancer per replicate, the 2,000 statistics sorted, and the bounds read at indices `int(0.025 * 2000)` and `int(0.975 * 2000)`.
* **Other seeds**: the sampling of the gene-recognition audit (`scripts/analysis/gene_recognition_audit/build_set.py`) and of the vocabulary readings uses seed 42; the proximity script names its own seed.

An independent recomputation that did not import these scripts reproduced the OncoDigger intervals (88.9-97.4% and 74.2-83.2%), the pooled-control gap (20.0 to 32.1 points) and the interaction intervals (21.1 to 41.6 and 31.1 to 47.9 points) once these draws were used. With a different seed or a per-replicate loop the bounds moved by up to about one point while every point estimate and p-value was unchanged.

## Known behaviours

* `.pptx` files are rewritten by the figure scripts with different zip-container bytes but identical slide content (checked part by part); they are therefore not hashed in `outputs_record.csv`.
* `retrieval_vocab_2x2_faithful_timing.csv` and `true_raw_flat_count_2x2_timing.csv` hold wall-clock seconds and differ on every run.
* `check_citations.py` prints a WARNING listing nine in-text citations and three references that its pattern cannot match: GeneReviews chapters cited by book date (Idos, Kamihara), "Research Network" for the TCGA Research Network, month names and the NLM name read as surnames ("September 2026", "Medicine 2026"), and four works whose in-text form the pattern does not parse (Kuehl, Metzler, Tao, Tuason). Each of them has an entry in the reference list; this was checked by hand.
* The evidence-audit worksheet (308 abstracts) was drawn from the 2026-09-27 ranking as it stood before the exact-background-count correction; its pair list is committed in `data/derived/evidence_audit/audited_top10_pairs.csv`.

## What was tested

Run on 2026-10-18 and 2026-10-19 from fresh clones, with the 2026-09-27 corpora:
* stages 00 to 40 (clone before the stage-order correction): every ranking file equal to the committed file; the only differences were four frequency-baseline files that had been stale from stage order (since corrected and no longer part of this release; no reported number used them);
* stages 45, 50, 70 (clone after the correction): pass; stage 70 reported that a second pass of stage 50 changed nothing;
* stage 65 (a later clone): all eight steps ran, exit 0; the two audits that need raw abstracts were skipped; no committed data file changed except a timing file;
* stage 60 on the final state: pass, including the number, citation and provenance checks.

Stages 00 to 40 were not repeated after the correction; the later changes touched only the audit script, the citation check, the provenance record and committed summaries.
