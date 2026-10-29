# OncoDigger: analysis code and data for the IntOGen comparison

Code, derived outputs and a build-provenance record for the manuscript

> OncoDigger Prioritizes Mutation-Associated Genes in Cancer-Specific PubMed Abstracts

which evaluates OncoDigger's mutation-focused gene rankings against IntOGen (Martínez-Jiménez et al., 2020)
as an external, cohort-computed reference. This repository is the code and data release for that manuscript.

**All results are computed from one corpus build, the 2026-09-27 release (1,640,505 abstracts).** See
`data/provenance/README.md`. A script that reads the corpora refuses to run on any other build.

## Reproducing the results

`scripts/reproduce/` holds one script per stage; running them in order, with the corpora, regenerates the rankings, tables,
statistics, figures and PDF of the manuscript, except the parts listed under "Not rerun" below. The full description (inputs, commands, measured times, what is not rerun, what was tested) is
in **`docs/REPRODUCTION.md`**.

```
export ONCODIGGER_CORPORA=/path/to/pubmed_processed     # the 19 databases of the 2026-09-27 build (about 27 GB)
bash scripts/reproduce/run_all.sh
```

What each part needs:

* **Stages 50, 60, 70 (tables, statistics, figures, PDF, fixed-point check)** need only this repository and the packages
  listed under "Install" (stage 60 rebuilds the figure sources with `python-pptx` and `pillow` and the PDF with `reportlab`).
* **Stage 45 (comparators)** also needs the DISEASES and PubTator cancer subsets (`data/external/diseases/`,
  `data/external/pubtator/`). They are in the Zenodo archive of this release, which is the complete package, and are not in the GitHub
  repository (the 257 MB DISEASES file exceeds GitHub's 100 MB file limit); copy them from the Zenodo archive into those folders. Stage 45 checks them against
  `data/external/PROVENANCE.md` and refuses any other file;
  CancerMine files are fetched from the fixed Zenodo record and checksum-verified.
* **Stages 00 to 40 and 65 (rankings, controls, the metastatic-literature exclusion test, corpus audits)** need the 2026-09-27 corpus
  databases, identified by SHA-256 in `data/provenance/corpus_build_2026-09-27.json`. They are not bundled in the
  repository. Two audit scripts additionally need the case-preserved raw abstracts and are skipped without them. The vocabulary screen in
  stage 65 fetches a public 20,000-word list when it is absent and checks it against the recorded SHA-256.
* **Not rerun by the stages:** the four scripts that query live services (their outputs are committed and
  checksummed), the PowerPoint-to-PNG figure export, and the gene-recognition audit (Supplementary Results S16, Supplementary
  Table S13). For the audit, the reader's identity verdicts are stored data that cannot be recomputed, the reading sets are not
  distributed, and the counterfactual harness needs the corpora; `data/derived/gene_recognition_audit/README.md` says which steps
  can be rerun (`analyze_cf.py` can, from the committed `cf_*.json` files).

## Contents

- `manuscript/` — the manuscript source (`OncoDigger_manuscript1.1.md`) and the rendered PDF.
- `scripts/reproduce/` — the numbered reproduction stages and `env.sh`; `scripts/analysis/`, `scripts/figures/`,
  `scripts/manuscript/` — analysis, figure and manuscript-building code; `scripts/analysis/gene_recognition_audit/` — the
  as-run scripts of the gene-recognition audit.
- `data/external/` — IntOGen compendium (release 2024.09.20, CC0), the intogen-plus whitelist (from the intogen-plus
  tag `v2024`), the saved ChatGPT response, and `PROVENANCE.md` for the external comparator datasets. The DISEASES and PubTator
  cancer subsets (`diseases/`, `pubtator/`) are in the Zenodo archive only; the CancerMine files are fetched by stage 45.
- `data/derived/`, `data/canonical/` — every derived table the manuscript cites, including `data/derived/gene_recognition_audit/`
  (per-credit evidence records with PubMed identifiers and titles, summaries and counterfactual runs). The canonical ranking is the
  main-pipeline ranking with the MET/REST/FH correction applied (`apply_vocab_correction_to_ranking.py`); the
  uncorrected version is kept as `precision_at_k_detail_pre_vocab_fix.csv.bak`.
- `data/supplementary_tables/` — the Repository Tables RT1 to RT12 cited in the manuscript, with `INDEX.md`.
- `data/provenance/` — corpus-build manifest, PMID manifests of the collections, raw-reconstruction record and the outputs record (see its README).
- `figures/source/`, `figures/png/` — editable PowerPoint sources and PNGs (2010 px wide, 300 DPI at the 6.7 in print width).
- `reference_engine_code/` — vendored engine source (commit `4e39fc0`, 2026-05-31) and vocabulary files (the KEGG
  and phytochemical lists are header-only placeholders, see its README), plus `gene_quality_test_engine/` (live-application commit `d592370`).
- `docs/` — `REPRODUCTION.md`, the exact-background specification, the vocabulary-policy history and two revision notes.
- `exploratory/background_rebuild/` — `hook/`, the exact-background count that every ranking in the paper uses (set by
  `BGFIX_ENABLE=1`), and the pilot that tested the matcher; its `README_superseded.md` says what is final and what is superseded.
  `exploratory/formula_study/` — the post hoc review of the relevance formula in which the background-count defect was found (Methods 3.2).

## Checking the repository

```
python -X utf8 scripts/analysis/check_provenance.py          # outputs unchanged, outputs agree with each other
ONCODIGGER_ROOT=<corpus folder> python -X utf8 scripts/analysis/check_provenance.py   # also verifies the corpus build
```

## Install

The tests recorded below used Python 3.12 and the versions in `requirements-used.txt` (without `comtypes`, which only the
PowerPoint-to-PNG export needs):

```
pip install pandas==3.0.6 numpy==2.2.6 scipy==1.14.1 requests==2.34.2 python-pptx==1.0.2 pillow==12.3.0 reportlab==4.5.1
```

`requirements.txt` holds an older set of pins; its comments explain why that set is installed in two steps.

## Verification status

<!-- VERIFICATION-STATUS: filled in only from runs that were actually executed -->
Reproduction tests of 2026-10-18 and 2026-10-19 were run from fresh clones of the development repository (its history is not part of
this release), with the 2026-09-27 corpora. The record is in `docs/REPRODUCTION.md` ("What was tested"):

* **Stages 00 to 40:** every ranking and control output equals the committed file. The only differences were four
  frequency-baseline files that had been stale because of stage order; they are corrected and no longer part of this release,
  and no reported number used them.
* **Stages 45, 50, 70:** pass; stage 70 shows that a second pass of stage 50 changes nothing.
* **Stage 65:** all eight audit steps ran and reproduce their committed outputs (a timing file differs); the two
  audits that need raw abstracts were skipped.
* **Stage 60:** figures, PDF and the number, citation and provenance checks pass.
* Stages 00 to 40 were not repeated after the stage-order correction; later changes touched only the evidence-audit
  script, the citation check, the provenance record and committed summaries.
* The evidence worksheet (`data/derived/evidence_audit/worksheet.csv`) is an unreviewed sample; its verdict
  columns are empty and the manuscript makes no claim about it.
* **This release as exported (2026-10-29)**, tested without the corpora, with Python 3.12 and the versions in `requirements-used.txt`: stages 45, 50, 70 and 60 pass on the unpacked archive (all files, including the DISEASES and PubTator
  subsets) and stages 50 and 70 pass on a clone of the GitHub repository's tracked files; stage 45 on that clone stops with a message that the two subsets are missing, as
  described above. Stage 70 in this release no longer depends on git. Nothing the stages write differs from the release files except the PDF and the `.pptx` zip containers
  (every part inside each `.pptx` is identical) and the CancerMine files fetched by stage 45. Stages 00 to 40 and 65 were not rerun on the exported release.

## License and citation

Code and data are released under the MIT License (see `LICENSE`). See `CITATION.cff` for how to cite this work.
