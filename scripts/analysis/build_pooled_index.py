"""
build_pooled_index.py — one lexical index over all 19 "Common cancers" corpora (no corpus separation), used for
the undifferentiated global control (all_pubmed_proxy.py; Methods 2.9).

Documents are deduplicated by PMID (a paper indexed in several cancer corpora appears once; its postings are
identical in each source because they derive from the same title/abstract). IDF and average document length are
recomputed over the pooled collection with the engine's formulas (scripts/build_lexical_index.py in the
OncoDigger repository): idf = ln(1 + (N - df + 0.5) / (df + 0.5)), avg_doc_length = mean(doc_length).

Output (outside this repository): <OncoDigger>/data/pubmed_processed/pooled_common_cancers_all/lexical_index.db
Run:    python -X utf8 scripts/analysis/build_pooled_index.py
"""
from __future__ import annotations

import math
import sqlite3
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_query_planes as rq  # noqa: E402

BASE = rq.DEFAULT_ONCODIGGER / "data" / "pubmed_processed"
OUT = BASE / "pooled_common_cancers_all" / "lexical_index.db"
SCHEMA = """
CREATE TABLE documents (pmid TEXT PRIMARY KEY, title TEXT, abstract TEXT, journal TEXT, journal_iso TEXT,
    publication_year INTEGER, doc_length INTEGER, mesh_terms TEXT, publication_types TEXT, authors TEXT,
    article_ids TEXT, publication_date TEXT);
CREATE TABLE postings (token TEXT NOT NULL, pmid TEXT NOT NULL, tf INTEGER NOT NULL);
CREATE TABLE idf (token TEXT PRIMARY KEY, idf REAL NOT NULL);
CREATE TABLE corpus_stats (avg_doc_length REAL NOT NULL, doc_count INTEGER NOT NULL);
"""


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    if OUT.exists():
        OUT.unlink()
    con = sqlite3.connect(str(OUT))
    con.executescript("PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF; PRAGMA temp_store=MEMORY;" + SCHEMA)
    t0 = time.time()
    for key, label in rq.CANCER_MAP:
        src = BASE / f"{key}_all" / "lexical_index.db"
        con.execute("ATTACH DATABASE ? AS src", (str(src),))
        con.execute("DROP TABLE IF EXISTS new_pmids")
        con.execute("CREATE TEMP TABLE new_pmids AS SELECT pmid FROM src.documents "
                    "WHERE pmid NOT IN (SELECT pmid FROM main.documents)")
        con.execute("CREATE INDEX temp.ix_new ON new_pmids(pmid)")
        con.execute("INSERT INTO main.documents SELECT d.* FROM src.documents d JOIN new_pmids n ON d.pmid = n.pmid")
        con.execute("INSERT INTO main.postings SELECT p.token, p.pmid, p.tf FROM src.postings p "
                    "JOIN new_pmids n ON p.pmid = n.pmid")
        n_new = con.execute("SELECT COUNT(*) FROM new_pmids").fetchone()[0]
        con.commit()
        con.execute("DETACH DATABASE src")
        print(f"{label:20} +{n_new:,} new documents  [{time.time() - t0:.0f}s]", flush=True)
    n = con.execute("SELECT COUNT(*) FROM documents").fetchone()[0]
    avgdl = con.execute("SELECT AVG(doc_length) FROM documents").fetchone()[0]
    con.execute("CREATE INDEX idx_postings_token ON postings(token)")
    con.executemany("INSERT INTO idf VALUES (?, ?)",
                    ((tok, math.log(1 + (n - df + 0.5) / (df + 0.5)))
                     for tok, df in con.execute("SELECT token, COUNT(*) FROM postings GROUP BY token").fetchall()))
    con.execute("INSERT INTO corpus_stats VALUES (?, ?)", (avgdl, n))
    con.commit()
    con.close()
    print(f"Pooled index: {n:,} documents, avg length {avgdl:.1f}  [{time.time() - t0:.0f}s]  -> {OUT}")


if __name__ == "__main__":
    main()
