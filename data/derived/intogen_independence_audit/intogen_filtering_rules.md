# IntOGen driver-selection rules relevant to benchmark independence

Sources (verbatim copies/extracts in `source_docs/`): IntOGen 2024 documentation (intogen-plus.readthedocs.io, pages `postprocessing`, `drivers_combination`, `methods_driver_discovery`, `data`, fetched 2026-10-11) and the `bbglab/intogen-plus` repository at tag **v2024**, commit `29e9ec5f6e5d1db1ffe39a0a2f3cab3d2ec4ee22` (2024-09-19), read but never executed. Released files audited: `Compendium_Cancer_Genes.tsv` and `Unfiltered_drivers.tsv` of release 2024.09.20 (folder stamp 2024-06-18; its README header still says "RELEASE 20230531").

**Version caveat.** The code tag post-dates the data-folder stamp. That the tag's logic reproduces the released tables is *tested* below (§ Consistency tests), not assumed.

Legend for "uses": **MUT** tumour mutation data; **CGC** Cancer Gene Census membership/annotation; **LIT** CancerMine literature; **OTHER** other external knowledge.

| # | Rule | What it does | Where (code) | Uses | Effect class | Can it affect an OncoDigger-positive gene's benchmark status? |
|---|---|---|---|---|---|---|
| 1 | Driver-detection methods (dNdScv, OncodriveFML, OncodriveCLUSTL, cBaSE, MutPanning, HotMAPS3D, smRegions) | per-cohort p-values from mutation data | `docs: methods_driver_discovery` | MUT | statistical detection | Source of the mutation evidence |
| 2 | Method weights for combination | weights per cohort chosen to maximise enrichment of CGC genes near the top of each method's ranking (CGC-Score), with a negative set of non-driver genes | `combination/.../grid_optimizer.py`, `schulze.py`; docs `drivers_combination` | MUT + **CGC** | knowledge-informed calibration of the statistic | Yes, for all genes, in a way that cannot be removed or isolated without rerunning the pipeline |
| 3 | Weighted Stouffer combination, q-value over all genes | `QVALUE_stouffer_w` = released `QVALUE_COMBINATION` | `stouffer_script.py::combine_pvals` | MUT (weights from #2) | statistic | Yes: basis of Tier 1/3 |
| 4 | CGC-only FDR | for genes in CGC the same combined p-values are FDR-corrected **among CGC genes only** (`QVALUE_CGC_stouffer_w`) | `stouffer_script.py` lines ~120-128 | MUT + **CGC** | structurally lower significance bar for CGC genes | Yes (Tier 2 rescue) |
| 5 | Tiers (pipeline uses the `.05` output: `intogen.nf` line ~586) | Tier 1 = ranked genes before the first gene with q>0.05; Tier 3 = q≤0.05 but beyond that block; Tier 4 = the rest; **Tier 2 = CGC rescue**: genes not in Tier 1 whose CGC-only q<0.25 | `create_tiers_drivers.py::classify_genes_tiers, rescue_genes, run` | MUT; Tier 2 uses **CGC** | classification / admission | Tier 2 rows exist only for CGC genes. |
| 6 | Candidate (`DRIVER`) rule | `DRIVER = (TIER≤3 and CGC gene)` **or** `(≥2 significant methods and TIER≤3)`; "significant method" = q<0.1 for that method | `vetting.py::get_drivers`; `stouffer_script.py::select_significant_bidders` | MUT; waiver uses **CGC** | candidate admission | A CGC gene needs no significant individual method; a non-CGC gene needs ≥2. |
| 7 | Filters applied to candidates, in order | (a) not expressed (TCGA, ≥80% negative); (b) hypermutation Signature 9 in lymphoid types; (c) **>2 mutations in one sample and gene not CGC**; (d) **germline-prone (MUTS/SAMPLE>1) and gene not CGC**; (e) olfactory receptor; (f) known artifact / suspected artifact lists; (g) **no CancerMine references, gene not CGC, gene not in the 80-gene whitelist → "Lack of literature evidence"**; else PASS | `vetting.py::perform_vetting` | MUT; (c)(d)(g) use **CGC**; (g) uses **LIT** (`n_papers` = count of CancerMine sentences for the gene, snapshot 2021-12-07, via `discovery.py::include_literature`) | final filtering (can only REMOVE candidates) | (c)(d)(g) exempt CGC genes; (g) can remove a literature-poor non-CGC gene; none of the filters can ADD a gene |
| 8 | CGC gene definition | `CGC_GENE` = gene present in the parsed CGC file. The parser's **Tier-1-only filter is commented out** (`build/datasets/cgc/parse.py` lines 114-115), so the file is not restricted to CGC Tier 1 | `build/datasets/cgc/parse.py`; `vetting.py::vet` | **CGC** | defines exemptions in #6 and #7 | Documentation (`postprocessing`) says "Tier1 CGC"; code does not enforce Tier 1 |
| 9 | CGC tumour-type annotation (3 levels) | annotates whether the gene is a CGC gene for that tumour type | `vetting.py::get_cancer_genes`; docs "Classification according to MSKCC oncotree" | CGC | annotation only | No (not used for admission) |
| 10 | Mode of action | inferred role combined with CGC role | `role.py`; docs | MUT + CGC | annotation only (`ROLE`) | No |
| 11 | Known long genes (TTN, OBSCN, RYR2...), olfactory receptors, artifact lists, black list | removes likely false positives | `filters.py`, `artifacts.json`, `black_listed.txt` | OTHER (prior knowledge) | final filtering | Only removes |
| 12 | `NUM_PAPERS` field | CancerMine sentence count per **gene** (identical across cancers for all 7,679 genes in the table) | `discovery.py` | LIT | annotation + input to #7(g) | Not cancer-specific |

## Direct answers to the audit questions

- **Is CancerMine literature required, optional, annotation-only or conditional?** *Conditionally used as an exclusion gate*: required (n_papers>0) only for genes that are not CGC and not whitelisted; irrelevant for CGC genes. It is also reported as the annotation `NUM_PAPERS`. It never admits a gene and never changes a q-value or tier.
- **Statistical detection vs knowledge-informed catalog construction.** Detection (#1) and the global q (#3) are mutation-based, but the method weights (#2) are CGC-calibrated. Catalog construction adds CGC-specific waiver (#6), CGC-only FDR rescue (#4/#5), CGC filter exemptions (#7c,d,g) and a CancerMine exclusion gate (#7g).
- **Counterfactual "remove only the literature criterion, hold everything else fixed."** Computable exactly from the released table: it ADDS the 142 rows currently labelled "Lack of literature evidence" (131 genes, all with ≥2 significant methods and zero CancerMine papers) and removes none. A counterfactual without CGC rules or weights requires re-executing the full pipeline on the raw cohorts and was **not** performed.

## Consistency tests of the v2024 logic against the released 2024.09.20 table (`verification.json`)

- TIER ∈ {1,3} if and only if `QVALUE_COMBINATION` ≤ 0.05: 0 exceptions in 23,444 rows → threshold 0.05 confirmed.
- All 193 catalog rows with q>0.05 are Tier 2 (CGC-rescue only).
- All 142 "Lack of literature", 168 ">2 mutations" and 21 "Germline" rows have ≥2 significant methods and none is Tier 2: exactly what the code predicts if those filters hit only non-CGC genes.
- Catalog rows (FILTER==PASS) = Compendium rows (4,478/4,478).
- No catalog row is Tier 4, expression-, artifact- or olfactory-flagged.

## What CGC membership we can and cannot establish

CGC membership is **not** in the released tables and the CGC file is licence-restricted. A catalog row is *certainly* a CGC gene (rule logic) when it passes with fewer than 2 significant methods, is Tier 2, passes the >2-mutations or germline filter while the flag condition holds, or passes with zero CancerMine papers while not whitelisted. This gives a **lower bound** (473 of 633 genes). All other genes are "CGC status unknown", never "non-CGC".
