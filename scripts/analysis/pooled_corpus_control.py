"""
pooled_corpus_control.py — undifferentiated global control for corpus separation (Methods 2.9; Figure 3C,
Supplementary Figure 6).

The identical OncoDigger pipeline (query 'mutation mutations', BM25 + enrichment, same settings, gene exclusions,
longest-match correction and hash seed as run_query_planes.py) is run once on a single index pooling all 19
cancer corpora (build_pooled_index.py). The result is one global gene ranking with no cancer-type separation; its
top-10 is compared with each cancer's own top-10 in build_primary_tables.py.

Output: data/derived/pooled_control_mutation.csv  (rank, gene), ranks 1-50
Run:    python -X utf8 scripts/analysis/pooled_corpus_control.py
"""
from __future__ import annotations

import csv
import os
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_query_planes as rq  # noqa: E402

POOLED = rq.DEFAULT_ONCODIGGER / "data" / "pubmed_processed" / "pooled_common_cancers_all" / "lexical_index.db"
OUT = rq.ROOT / "data" / "derived" / "pooled_control_mutation.csv"


def main() -> None:
    if os.environ.get("PYTHONHASHSEED") != "0":
        import subprocess
        env = dict(os.environ, PYTHONHASHSEED="0")
        sys.exit(subprocess.run([sys.executable, "-X", "utf8"] + sys.argv, env=env).returncode)
    root = rq.DEFAULT_ONCODIGGER
    sys.path.insert(0, str(root / "src"))
    from oncodigger.search import discovery
    from oncodigger.search.discovery import associated_entity_tables
    from oncodigger.search.lexical import bm25_search, normalize_query
    from oncodigger.search.lexicons import load_lexicons
    genes = rq.longest_match_lookups(rq.suppress_ambiguous_symbols(load_lexicons(root)["Genes"]), discovery)
    engine = (bm25_search, normalize_query, associated_entity_tables, {"Genes": genes})
    con = rq.open_corpus(POOLED, discovery)
    n = con.execute("SELECT doc_count FROM corpus_stats").fetchone()[0]
    ranked = rq.rank_genes(con, rq.PLANES["mutation"], engine)
    con.close()
    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["rank", "gene"])
        w.writerows(enumerate(ranked, 1))
    print(f"Pooled corpus: {n:,} documents; global top-10: {ranked[:10]}")


if __name__ == "__main__":
    main()
