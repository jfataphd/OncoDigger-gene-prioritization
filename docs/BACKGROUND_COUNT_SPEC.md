# Specification: background count in the gene enrichment term

Drafted 2026-10-15, before any exact recomputation. Status: APPROVED by the author 2026-10-15 (Stages 0 to 2 authorised; the full rerun needs a separate go-ahead after the pilot is reviewed). Nothing here changes the released paper, the
submitted repository or Zenodo.

## The defect

The relevance score is `R = W * ln(1 + E) * ln(1 + S)` with `E = (S / N) / (B / C)`, where S is the number of retrieved papers (pool of N)
in which the gene is found, B the number of papers in the cancer collection (C papers) in which the gene occurs, and W the rank-weighted sum
over the S papers. The manuscript describes B as the gene's background frequency in the cancer collection.

The engine (`association.py`, `entity_background_papers`) does not compute that quantity:

1. It takes the smallest postings count among the gene's matched single-word terms (the rarest alias that happened to occur in the pool),
   not the number of papers containing any alias.
2. Matched multi-word terms contribute nothing, so a gene matched only as a phrase (for example NKX2-1, normalised `nkx2 1`) has B = 0, and the
   caller then substitutes B = S. For such genes E reduces to C / N whatever the gene's real frequency.
3. Postings tokens keep hyphens and the pool matcher does not, so hyphenated forms are in S but not in B.
4. The retrieval postings are not the same recogniser as the pool matcher (no longest-match masking, no phrase matching, no alias suppression).

Evidence (exploratory, `exploratory/formula_study/`): the engine's B is smaller than a postings-based union over aliases for 1,553 of 3,892
candidate genes and more than twice as small for 688. Example: NKX2-1 in Lung Cancer is reported "found in 14 of 1000 relevant papers and
14 of 178,908 corpus papers".

## Definition (to be implemented)

For each cancer collection and each gene entity, `B` = the number of distinct PMIDs in the collection in which the SAME matcher that produces
S finds that entity. The matcher is the one in `associated_entity_tables`: normalised title + abstract text (longest-match masking applied),
single-word lookups on document tokens, phrase lookups by anchor token and substring test, the active alias suppression, and the redaction of
the MET, REST and FH false positives already present in the corrected database text. Counted once per PMID.

`S` and `B` therefore use one definition of "gene found in a paper", and `S <= B <= C` must hold for every candidate gene.

Entity level: the unit is the lexicon entity id, as for S. The pipeline's remap (P53 to TP53, HER2 to ERBB2) and its first-row deduplication are
left as they are, so S and B for a remapped name are those of the entity that is kept.

## What stays fixed

BM25 settings and multipliers, the pool of up to 1,000 papers, the query, the cleaned vocabulary and its corrections, W, and the formula
`R = W * ln(1 + E) * ln(1 + S)`. Only B changes. The enrichment floor `background = max(background, supporting)` stays as a guard but is
expected never to bind; the number of times it binds is reported.

Also needed: B for every other index the pipeline ranks on (the pooled 19-collection index, the 2022 and 2023 truncated collections, the
true-raw and uncleaned-vocabulary runs, which have their own recogniser), each computed with the matcher of that run.

## Tests (must pass before any rerun is used)

1. Replica check: counting the matcher over the retrieved pool reproduces the engine's S for every candidate gene in every cancer (exact).
2. Invariant: `S <= B <= C` for every candidate gene in every cancer, no exceptions.
3. Phrase-only genes (NKX2-1, HLA-A, HLA-B, H3-3A) get a B that counts phrase matches.
4. Merged ids (TP53/P53, ERBB2/HER2) are counted per entity as S is.
5. Hyphenated and unhyphenated forms (HER-2, HER2) both counted.
6. Determinism: two runs give identical counts (`PYTHONHASHSEED=0`).
7. For genes whose B is unchanged by the fix (engine B equals exact B) the relevance is unchanged to the printed precision.

## How the result is used

The corrected background becomes the primary analysis once the tests pass, whichever direction concordance moves. The previous results stay
available as git tag `pre-background-fix-2026-10-15`. Approximate figures from the exploratory study (144/190 exact with an alias-union
approximation) are not used. The relevance formula is not changed or tuned in this step. The way the defect was found is disclosed in the
manuscript.

## Refinements adopted at approval

1. Candidate universe. B is computed for every entity the matcher finds in the collection, not only the current top 50, so that a corrected
   ranking can promote genes from below rank 50 with a correct B.
2. One recogniser for S and B. Normalisation, disambiguation, longest-match masking, aliases, identity mapping and per-PMID deduplication are
   shared. If a PMID contributes to S it contributes to B. The matcher runs over the full vocabulary of competing entities, because
   longest-match masking depends on which other names are present. A gene is never matched in isolation.
3. Direct matching is canonical. Any inverted-index shortcut is allowed only if it reproduces the full-text matcher exactly. A small
   hand-built corpus with known expected counts is part of the test set.
4. Version-stamped cache. Every cached count is stored with the corpus build hash, vocabulary hash, matcher version and redaction policy, and is
   refused when any of them differs.
5. Dependency split. Analyses that use B are rerun. External rankings and the flat publication-count baselines do not use B and must come out
   unchanged. If they change, the cause is investigated before anything else is used.
6. Extra invariant. For every gene with S > 0, `0 < E <= C / N`, with N the number of distinct papers in the pool.

## Frozen decision (recorded before the pilot)

The original relevance formula is retained. The background-count defect is corrected whether agreement with IntOGen rises or falls. No
formula or parameter is chosen on the basis of the corrected benchmark results. The old results (142/190 exact, 177/190 any-cancer) are kept as
a historical comparison that documents the effect of the correction and are not alternative primary findings.

## Pilot (Stages 0 to 2) and its gate

Cancers: Myeloma (small, 26,591 papers) and Lung Cancer (large, 178,908 papers, prior MET work). The pilot answers: are the counts accurate
(tests 1 to 7), is the run time practical for 19 cancers plus the supplementary runs, and can it be reproduced independently. It does not decide
whether the correction is adopted. The full rerun is authorised separately after the pilot is reviewed.
