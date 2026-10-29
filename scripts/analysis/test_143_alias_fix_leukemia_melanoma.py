"""
test_143_alias_fix_leukemia_melanoma.py -- targeted check, NOT part of the frozen 2x2
decomposition. Tests whether gene_quality.py's 143-alias-collision removal (confirmed
NOT applied by this paper's own frozen pipeline, run_query_planes.py) changes Leukemia's
or Melanoma's "mutation, mutations" top-10, since those are the two cancers the user
recalled seeing a substantial Raw-vs-Vetted difference for on the live site.

Compares three vocabularies, same real BM25 retrieval pool and real relevance formula
throughout (only the vocabulary varies):
  A. This paper's current "cleaned" vocabulary (load_lexicons -> suppress_ambiguous_symbols
     -> longest_match_lookups) -- what every number in the manuscript is based on.
  B. gene_quality.py's "raw" output (no cleaning at all).
  C. gene_quality.py's "vetted" output (the 143-alias-collision removal PLUS the 20-symbol-
     style exclusions gene_quality.py implements under its own, independently-written policy).

Run: python -X utf8 scripts/analysis/test_143_alias_fix_leukemia_melanoma.py
"""
from __future__ import annotations

import sys
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ONCODIGGER = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
# Uses a separate, consistently-versioned copy of the engine (commit d592370, the one
# gene_quality.py actually belongs to -- the live site's current deployed commit, per
# user confirmation it hasn't been updated since). Deliberately NOT reference_engine_code/,
# which stays pinned to commit 4e39fc06 for the manuscript's own reproducibility claim.
sys.path.insert(0, os.environ.get("GENE_QUALITY_TEST_ENGINE", str(ROOT / "reference_engine_code" / "gene_quality_test_engine")))

import sqlite3

from oncodigger_search import gene_quality as gq
from oncodigger_search import discovery
from oncodigger_search.discovery import associated_entity_tables
from oncodigger_search.lexical import bm25_search, normalize_query
from oncodigger_search.lexicons import load_lexicons

QUERY = "mutation, mutations"
SEARCH_POOL, ANALYSIS_PAPERS = 1250, 1000
OPTIONS = {"year_min": 2000, "year_max": 2026, "title_match_multiplier": 1.2,
           "mesh_match_multiplier": 1.15, "recency_boost": True, "discovery_mode": "Primary Evidence"}
GENE_EXCLUDE = {"MICE", "DCR", "PC", "C2", "FBN1", "SGCG", "WDHD1"}
SUPPRESS_SYMBOLS = {"GC", "HR", "HCCS", "BPIFA4P", "MB", "ELN", "CP", "UBC", "CAMP", "EFS", "IVD", "SCT", "SON"}
GENE_REMAP = {"P53": "TP53", "HER2": "ERBB2"}

CANCERS = [("leukemia_cancer_inclusive", "Leukemia"), ("melanoma_cancer_inclusive", "Melanoma")]


def suppress_ambiguous_symbols(genes_lexicon: dict) -> dict:
    single = {k: v for k, v in genes_lexicon["single"].items() if k.upper() not in SUPPRESS_SYMBOLS}
    return {"single": single, "phrases": genes_lexicon["phrases"]}


def rank_genes(con, vocab_lexicon: dict) -> list[str]:
    discovery._sqlite_doc_counts.clear()
    discovery._sqlite_caches.clear()
    tokens = normalize_query(QUERY)
    hits = bm25_search(con, tokens, SEARCH_POOL, OPTIONS)[:ANALYSIS_PAPERS]
    tables = associated_entity_tables(con, hits, QUERY, (2000, 2026), tokens, {"Genes": vocab_lexicon},
                                      limit=100, discovery_mode="Primary Evidence", corpus_query="")
    genes_df = tables.get("Genes")
    ranked, seen = [], set()
    if genes_df is None or genes_df.empty:
        return ranked
    for _, row in genes_df.iterrows():
        gene = GENE_REMAP.get(row["Gene"], row["Gene"])
        if gene in GENE_EXCLUDE or gene in seen:
            continue
        seen.add(gene)
        ranked.append(gene)
        if len(ranked) >= 10:
            break
    return ranked


def main() -> None:
    print("Building vocabularies...")
    lex_manuscript_raw = load_lexicons(DEFAULT_ONCODIGGER)["Genes"]
    lex_manuscript_clean = suppress_ambiguous_symbols(lex_manuscript_raw)

    raw_gq, vetted_gq, audit = gq.gene_lookups(DEFAULT_ONCODIGGER)
    print(f"gene_quality.py vetted audit: {len(audit)} logged actions "
          f"({sum(1 for a in audit if a['Action']=='Removed alias')} alias collisions removed)")

    for corpus_key, cancer in CANCERS:
        db = DEFAULT_ONCODIGGER / "data" / "pubmed_processed" / f"{corpus_key}_all" / "lexical_index.db"
        con = sqlite3.connect(str(db))
        print(f"\n=== {cancer} ===")
        top_a = rank_genes(con, lex_manuscript_clean)
        print(f"  A. This paper's cleaned vocabulary:  {top_a}")
        top_b = rank_genes(con, raw_gq)
        print(f"  B. gene_quality.py RAW:              {top_b}")
        top_c = rank_genes(con, vetted_gq)
        print(f"  C. gene_quality.py VETTED (143-fix):  {top_c}")
        print(f"  A vs C differ: {set(top_a) != set(top_c)}  -- diff: {set(top_a) ^ set(top_c)}")
        print(f"  B vs C differ: {set(top_b) != set(top_c)}  -- diff: {set(top_b) ^ set(top_c)}")


if __name__ == "__main__":
    main()
