"""
run_query_planes.py — cleaned OncoDigger rankings for any query across the 19 corpora.

Reproduces the post-alias-cleanup pipeline used for the primary mutation ranking
(data/canonical/precision_at_k_detail.csv) and applies it, with identical settings, to
the additional query planes (Methods 2.10) and the role-specific queries (Section 3.5):

  mutation            "mutation mutations"            (validation: must reproduce the main ranking)
  fusion              "fusion translocation"
  methylation         "methylation epigenetic"
  amplification       "amplification overexpression"
  tumor_suppressor    "tumor suppressor"
  oncogene            "oncogene"

Settings (identical to the primary analysis): Primary Evidence mode, years 2000-2026, recency
boost, title x1.2 / MeSH x1.15 multipliers, review down-ranking, top 1,000 of 1,250 BM25 matches,
alias-cleaned HGNC 2026 vocabulary, P53->TP53 and HER2->ERBB2, exclusion list below, and genes
with no corpus postings skipped.

NOTE: the public site (oncodigger.com) defaults to "Overview" mode rather than this script's
(and the manuscript's) "Primary Evidence" mode, and since 2026-09-27 shows both "Raw" and
"Vetted" gene results side by side rather than serving only the uncleaned engine as earlier
versions of this comment said; results from this script (and the manuscript) will still differ
from the live site's default view.

Requires the local OncoDigger code and SQLite corpora (not in this repository):
  --oncodigger-root  (default: the author's local OncoDigger working copy)

Output: data/derived/query_planes/<plane>.csv  (cancer, rank, gene), ranks 1-50

Corpora: the 19 "Common cancers" corpora (*_inclusive; 2000-2026, reviews included and down-weighted).
Python hash seed is fixed at 0 (the script re-executes itself) so near-tied genes rank identically on rerun.

Run:  python -X utf8 scripts/analysis/run_query_planes.py [--planes mutation fusion ...]
"""
from __future__ import annotations

import argparse
import os
import csv
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ONCODIGGER = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
OUT_DIR = ROOT / "data" / "derived" / "query_planes"

PLANES = {
    "mutation": "mutation mutations",
    "fusion": "fusion translocation",
    "methylation": "methylation epigenetic",
    "amplification": "amplification overexpression",
    "tumor_suppressor": "tumor suppressor",
    "oncogene": "oncogene",
    # Reviewer round 2, item 1e (docs/reviewer_revision_plan.md): secondary sensitivity
    # analysis only, never the primary benchmark -- using "driver" in the query risks
    # circularity with the driver-curation concept of the reference catalogues, the same concern already raised
    # in this paper about querying "tumor suppressor"/"oncogene" directly.
    "driver_mutation": "driver mutation driver mutations",
}
TOP_K, LIMIT, ANALYSIS_PAPERS, SEARCH_POOL = 50, 100, 1000, 1250
# The exclusion and suppression lists live in vocab_policy.py (one frozen, versioned policy for every analysis that uses the cleaned
# vocabulary); each module gets its own mutable copy so that robustness tests can lift an entry without touching the policy.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import vocab_policy as _vp  # noqa: E402

GENE_EXCLUDE = set(_vp.GENE_EXCLUDE)            # 7 genes with no safe alias: dropped entirely (Methods 3.3)
SUPPRESS_SYMBOLS = set(_vp.SUPPRESS_SYMBOLS)    # 13 bare symbols (policy v3) + SHE, RHO, CCS, NODAL (policy v4): term-level suppression
SUPPRESS_ALIASES = {k: set(v) for k, v in _vp.SUPPRESS_ALIASES.items()}   # RANK (TNFRSF11A), nudE (NDE1), Lobe (AKT1S1)
GENE_REMAP = {"P53": "TP53", "HER2": "ERBB2"}
CANCER_MAP = [
    ("breast_cancer_inclusive", "Breast Cancer"), ("lung_cancer_inclusive", "Lung Cancer"),
    ("colorectal_cancer_inclusive", "Colorectal Cancer"), ("prostate_cancer_inclusive", "Prostate Cancer"),
    ("melanoma_cancer_inclusive", "Melanoma"), ("bladder_cancer_inclusive", "Bladder Cancer"),
    ("kidney_cancer_inclusive", "Kidney Cancer"), ("pancreas_cancer_inclusive", "Pancreatic Cancer"),
    ("liver_cancer_inclusive", "Liver Cancer"), ("stomach_cancer_inclusive", "Stomach Cancer"),
    ("esophageal_cancer_inclusive", "Esophageal Cancer"), ("ovarian_cancer_inclusive", "Ovarian Cancer"),
    ("endometrial_cancer_inclusive", "Endometrial Cancer"), ("cervical_cancer_inclusive", "Cervical Cancer"),
    ("thyroid_cancer_inclusive", "Thyroid Cancer"), ("brain_cancer_inclusive", "Brain Cancer"),
    ("leukemia_cancer_inclusive", "Leukemia"), ("lymphoma_cancer_inclusive", "Lymphoma"),
    ("myeloma_cancer_inclusive", "Myeloma"),
]
OPTIONS = {
    "year_min": 2000, "year_max": 2026, "exclude_reviews": False, "downrank_reviews": True,
    "title_match_multiplier": 1.2, "mesh_match_multiplier": 1.15, "recency_boost": True,
    "discovery_mode": "Primary Evidence",
}


def longest_match_lookups(genes_lexicon: dict, discovery) -> dict:
    """Longest-match correction for nested gene full names (manuscript pipeline only).

    The engine credits a gene whenever its full-name phrase occurs in the text, with no longest-match rule, so a
    name nested inside another gene's longer name is matched inside that name (e.g. EGF's "epidermal growth
    factor" inside EGFR's "epidermal growth factor receptor"). Here every longer name that contains another
    gene's name is masked in the document text (spaces -> underscores, which the tokenizer still splits), and the
    masked form is registered as a phrase for the longer name's own gene(s). Symbols and aliases are unaffected.
    Returns the patched lookup; patches discovery.full_document_text in place.
    """
    import re
    names: dict[str, list] = {}
    for anchor, lst in genes_lexicon["phrases"].items():
        for term, entity in lst:
            names.setdefault(term, []).append(entity)
    ids = {t: {e["id"] for e in es} for t, es in names.items()}
    containers = set()
    for q in names:
        w = q.split()
        for i in range(len(w)):
            for j in range(i + 1, len(w) + 1):
                sub = " ".join(w[i:j])
                if (i, j) != (0, len(w)) and sub in ids and ids[sub] - ids[q]:
                    containers.add(q)
    pattern = re.compile(r"(?<![a-z0-9'])(?:" + "|".join(re.escape(q) for q in sorted(containers, key=len, reverse=True))
                         + r")(?![a-z0-9'])")
    phrases = {k: list(v) for k, v in genes_lexicon["phrases"].items()}
    for q in containers:
        for e in names[q]:
            phrases.setdefault(q.split()[0], []).append((q.replace(" ", "_"), e))
    original = discovery.full_document_text

    def masked_text(doc: dict) -> str:
        source = " ".join([doc.get("title", ""), doc.get("abstract", "")])
        if doc.get("_manuscript_masked_source") != source:
            doc["_oncodigger_full_text"] = pattern.sub(lambda m: m.group(0).replace(" ", "_"), original(doc))
            doc["_manuscript_masked_source"] = source
        return doc["_oncodigger_full_text"]

    discovery.full_document_text = masked_text
    print(f"longest-match correction: {len(containers)} container names masked", flush=True)
    return {"single": genes_lexicon["single"], "phrases": phrases}


def suppress_ambiguous_symbols(genes_lexicon: dict) -> dict:
    """Term-level suppression of ambiguous gene terms (manuscript pipeline only; policy in vocab_policy.py).

    Mirrors the live app's gene_quality.py SAFE_ALIASES mechanism: for each symbol in SUPPRESS_SYMBOLS removes only the bare, case-normalized
    symbol as a matchable single-word term pointing to that gene's entity, and for each (gene, alias) in SUPPRESS_ALIASES removes only that
    alias term for that gene. Every other term of the same gene (full approved name, other aliases, multi-word phrases) is untouched and
    still scores that gene's identity. If a token is shared with a different gene's entity, only the matching gene's entity is removed.
    Four of the 13 original symbols (GC, HR, HCCS, BPIFA4P) are already indexed under a safer alias by lexicons.py (preferred_alias column),
    so for them this is a verified no-op.
    """
    return _vp.apply_suppression(genes_lexicon, SUPPRESS_SYMBOLS, SUPPRESS_ALIASES)


def open_corpus(db_path, discovery) -> sqlite3.Connection:
    """Open a corpus index and reset the engine's per-connection caches.

    The engine caches corpus size (_sqlite_doc_counts, used as N in the enrichment score) and term checks
    (_sqlite_caches) keyed by id(connection). A new connection can reuse the address of a closed one and would
    then inherit the previous corpus's size; clearing the caches on every open prevents this.
    """
    discovery._sqlite_doc_counts.clear()
    discovery._sqlite_caches.clear()
    return sqlite3.connect(str(db_path), check_same_thread=False)


def rank_genes(con, query_text, engine) -> list[str]:
    bm25_search, normalize_query, associated_entity_tables, gene_lookups = engine
    tokens = normalize_query(query_text)
    hits = bm25_search(con, tokens, SEARCH_POOL, OPTIONS)[:ANALYSIS_PAPERS]
    tables = associated_entity_tables(con, hits, query_text, (2000, 2026), tokens, gene_lookups,
                                      limit=LIMIT, discovery_mode="Primary Evidence", corpus_query="")
    genes_df = tables.get("Genes")
    ranked, seen = [], set()
    if genes_df is None or genes_df.empty:
        return ranked
    for _, row in genes_df.iterrows():
        gene = GENE_REMAP.get(row["Gene"], row["Gene"])
        if gene in GENE_EXCLUDE or gene in seen:
            continue
        try:
            n = con.execute("SELECT COUNT(DISTINCT pmid) FROM postings WHERE token = ?", (gene.lower(),)).fetchone()[0]
        except Exception:
            n = 0
        if n == 0:
            continue
        seen.add(gene)
        ranked.append(gene)
        if len(ranked) >= TOP_K:
            break
    return ranked


def main() -> None:
    # The engine sums scores while iterating over sets; with Python's randomised string
    # hashing, near-tied genes can swap order between runs. Fix the seed for reproducibility.
    if os.environ.get("PYTHONHASHSEED") != "0":
        import subprocess
        env = dict(os.environ, PYTHONHASHSEED="0")
        sys.exit(subprocess.run([sys.executable, "-X", "utf8"] + sys.argv, env=env).returncode)
    ap = argparse.ArgumentParser()
    ap.add_argument("--oncodigger-root", type=Path, default=DEFAULT_ONCODIGGER)
    ap.add_argument("--planes", nargs="+", default=list(PLANES), choices=list(PLANES))
    ap.add_argument("--cancers", nargs="+", help="subset of cancer labels (testing; writes to --out-dir)")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--no-longest-match", action="store_true", help="disable the nested-name correction")
    ap.add_argument("--allow-nonstandard-build", action="store_true",
                    help="skip the corpus-build check (only for deliberately different corpora, e.g. the year-truncated cutoff runs)")
    args = ap.parse_args()

    # Refuse to run on anything other than the recorded corpus build. An earlier version of this script
    # silently fell back to a different default corpus folder and produced rankings from another build.
    if not args.allow_nonstandard_build:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from provenance import assert_corpus_build
        assert_corpus_build(args.oncodigger_root)

    unknown = sorted(set(args.cancers or []) - {label for _, label in CANCER_MAP})
    if unknown:
        sys.exit(f"--cancers: unknown label(s) {unknown}; valid labels: {[label for _, label in CANCER_MAP]}")

    sys.path.insert(0, str(args.oncodigger_root / "src"))
    from oncodigger.search import discovery
    from oncodigger.search.discovery import associated_entity_tables
    from oncodigger.search.lexical import bm25_search, normalize_query
    from oncodigger.search.lexicons import load_lexicons
    lexicons = load_lexicons(args.oncodigger_root)
    genes_suppressed = suppress_ambiguous_symbols(lexicons["Genes"])
    genes = genes_suppressed if args.no_longest_match else longest_match_lookups(genes_suppressed, discovery)
    engine = (bm25_search, normalize_query, associated_entity_tables, {"Genes": genes})

    db_base = args.oncodigger_root / "data" / "pubmed_processed"
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    for plane in args.planes:
        rows = []
        for key, label in CANCER_MAP:
            if args.cancers and label not in args.cancers:
                continue
            db = next((p for p in (db_base / key / "lexical_index.db", db_base / f"{key}_all" / "lexical_index.db")
                       if p.exists()), None)
            if db is None:
                sys.exit(f"missing index for {label}")
            con = open_corpus(db, discovery)
            ranked = rank_genes(con, PLANES[plane], engine)
            con.close()
            rows += [{"cancer": label, "rank": i, "gene": g} for i, g in enumerate(ranked, 1)]
            print(f"  {plane:16} {label:20} {' '.join(ranked[:10])}", flush=True)
        out = out_dir / f"{plane}.csv"
        with out.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["cancer", "rank", "gene"])
            w.writeheader()
            w.writerows(rows)
        print(f"Saved {out} ({len(rows)} rows)", flush=True)


if __name__ == "__main__":
    main()
